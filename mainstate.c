/* Copying and distribution of this file, with or without modification,
 * are permitted in any medium without royalty. This file is offered as-is,
 * without any warranty.
 */

/*! @file mainstate.c
 * @brief Main State machine for template application.
 *
 * Makes use of Framework HSM module.
	************************************************************************/

#include "template.h"
#include "mainstate.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

const Msg mainStateMsg[] = {
	{ FRAMESEQ_EVT },
	{ FRAMEPAR_EVT },
	{ IPC_GET_APP_STATE_EVT },
	{ IPC_GET_THRESHOLD_IMG_EVT },
	{ IPC_GET_RAW_IMG_EVT },
	{ IPC_SET_SHOWIMG_MODE_EVT }
};

/*********************************************************************//*!
 * @brief Inline function to throw an event to be handled by the statemachine.
 *
 * @param pHsm Pointer to state machine
 * @param evt Event to be thrown.
 *//*********************************************************************/
void ThrowEvent(struct MainState *pHsm, unsigned int evt)
{
	const Msg *pMsg = &mainStateMsg[evt];
	HsmOnEvent((Hsm*)pHsm, pMsg);
}

/*********************************************************************//*!
 * @brief Checks for IPC events, schedules their handling and
 * acknowledges any executed ones.
 *
 * @param pMainState Initalized HSM main state variable.
 * @return 0 on success or an appropriate error code.
 *//*********************************************************************/
static OSC_ERR HandleIpcRequests(MainState *pMainState)
{
	OSC_ERR err;
	uint32 paramId;
	struct IPC_DATA *pIpc = &data.ipc;
	struct OSC_IPC_REQUEST *pReq = &pIpc->req;

	err = CheckIpcRequests(&paramId);
	if (err == SUCCESS)
	{
		/* We have a request. See to it that it is handled
		 * depending on the state we're in. */
		switch(paramId)
		{
		case GET_APP_STATE:
			/* Request for the current state of the application. */
			ThrowEvent(pMainState, IPC_GET_APP_STATE_EVT);
			break;
		case GET_COLOR_IMG:
			/* Request for the live image. */
			ThrowEvent(pMainState, IPC_GET_THRESHOLD_IMG_EVT);
			break;
		case GET_RAW_IMG:
			/* Request for the live image. */
			ThrowEvent(pMainState, IPC_GET_RAW_IMG_EVT);
			break;
		case SET_CAPTURE_MODE:
			/* Set the debayering option. */
			ThrowEvent(pMainState, IPC_SET_SHOWIMG_MODE_EVT);
			break;
		case SET_EXPOSURE_TIME:
			// a new exposure time was given
			if(data.ipc.state.nExposureTime != *((int*)pReq->pAddr))
			{
				data.nExposureTimeChanged = true;
				data.ipc.state.nExposureTime = *((int*)pReq->pAddr);
			}
			data.ipc.enReqState = REQ_STATE_ACK_PENDING;//we return immediately
			break;
		case SET_THRESHOLD:
			// a new exposure time was given
			if(data.ipc.state.nThreshold != *((int*)pReq->pAddr))
			{
				data.ipc.state.nThreshold = *((int*)pReq->pAddr);
			}
			data.ipc.enReqState = REQ_STATE_ACK_PENDING;//we return immediately
			break;
		default:
			OscLog(ERROR, "%s: Unkown IPC parameter ID (%d)!\n", __func__, paramId);
			data.ipc.enReqState = REQ_STATE_NACK_PENDING;
			break;
		}
	}
	else if (err == -ENO_MSG_AVAIL)
	{
		/* No new message available => do nothing. */
	}
	else
	{
		/* Error.*/
		OscLog(ERROR, "%s: IPC request error! (%d)\n", __func__, err);
		return err;
	}

	/* Try to acknowledge the new or any old unacknowledged
	 * requests. It may take several tries to succeed.*/
	err = AckIpcRequests();
	if (err != SUCCESS)
	{
		OscLog(ERROR, "%s: IPC acknowledge error! (%d)\n", __func__, err);
	}
	return err;
}

Msg const *MainState_top(MainState *me, Msg *msg)
{
	struct APPLICATION_STATE *pState;
	switch (msg->evt)
	{
	case START_EVT:
		STATE_START(me, &me->showRaw);
		data.nExposureTimeChanged = true;
		data.ipc.state.nExposureTime = 25;
		data.ipc.state.nStepCounter = 0;
		data.ipc.state.nThreshold = 30;
		return 0;
	case IPC_GET_APP_STATE_EVT:
		/* Fill in the response and schedule an acknowledge for the request. */
		pState = (struct APPLICATION_STATE*)data.ipc.req.pAddr;
		memcpy(pState, &data.ipc.state, sizeof(struct APPLICATION_STATE));

		data.ipc.enReqState = REQ_STATE_ACK_PENDING;
		return 0;
	case FRAMESEQ_EVT:
		/* Timestamp the capture of the image. */
		data.ipc.state.imageTimeStamp = OscSupCycGet();
		data.ipc.state.bNewImageReady = TRUE;
		/* Sleep here for a short while in order not to violate the vertical
		 * blank time of the camera sensor when triggering a new image
		 * right after receiving the old one. This can be removed if some
		 * heavy calculations are done here. */
		usleep(4000);
		return 0;
	case FRAMEPAR_EVT:
	{
		/* we have a new image increase counter: here and only here! */
		data.ipc.state.nStepCounter++;
		/* debayer the image first -> to half size*/
		OscVisDebayerGreyscaleHalfSize( data.pCurRawImg, OSC_CAM_MAX_IMAGE_WIDTH, OSC_CAM_MAX_IMAGE_HEIGHT, ROW_BGBG, data.u8TempImage[GRAYSCALE]);
		/* Process the image. */
		if(data.ipc.state.enAppMode == APP_CAPTURE_COLOR)
			/* the parameter is not really required */
			ProcessFrame(data.u8TempImage[GRAYSCALE]);

		return 0;
	}
	case IPC_GET_THRESHOLD_IMG_EVT:
	case IPC_GET_RAW_IMG_EVT:
	case IPC_SET_SHOWIMG_MODE_EVT:
		/* If the IPC event is not handled in the actual substate, a negative acknowledge is returned by default. */
		data.ipc.enReqState = REQ_STATE_NACK_PENDING;
	return 0;
	}
	return msg;
}

Msg const *MainState_ShowThreshold(MainState *me, Msg *msg)
{
	bool bShowThreshold;

	switch (msg->evt)
	{
	case ENTRY_EVT:
		data.ipc.state.enAppMode = APP_CAPTURE_COLOR;
		data.pCurRawImg = data.u8FrameBuffers[0];
		return 0;
	case IPC_GET_THRESHOLD_IMG_EVT:
	{
		/* Write out the image to the address space of the CGI. */
		memcpy(data.ipc.req.pAddr, data.u8TempImage[THRESHOLD], sizeof(data.u8TempImage[THRESHOLD]));

		data.ipc.state.bNewImageReady = FALSE;

		/* Mark the request as executed, so it will be acknowledged later. */
		data.ipc.enReqState = REQ_STATE_ACK_PENDING;
		return 0;
	}
	case IPC_SET_SHOWIMG_MODE_EVT:
		/* Read the option from the address space of the CGI. */
		bShowThreshold = *((bool*)data.ipc.req.pAddr);
		if (bShowThreshold == FALSE)
		{
			/* Need to capture raw images from now on, this is done in the showRaw state.  */
			STATE_TRAN(me, &me->showRaw);
		}
		data.ipc.enReqState = REQ_STATE_ACK_PENDING;
		return 0;
	}
	return msg;
}

Msg const *MainState_ShowRaw(MainState *me, Msg *msg)
{
	bool bShowThreshold;

	switch (msg->evt)
	{
	case ENTRY_EVT:
		data.ipc.state.enAppMode = APP_CAPTURE_RAW;
		data.pCurRawImg = data.u8FrameBuffers[0];
		return 0;
	case IPC_GET_RAW_IMG_EVT:
	{
		/* Write out the current gray image to the address space of the CGI. */
		memcpy(data.ipc.req.pAddr, data.u8TempImage[GRAYSCALE], sizeof(data.u8TempImage[GRAYSCALE]));

		data.ipc.state.bNewImageReady = FALSE;

		/* Mark the request as executed, so it will be acknowledged later. */
		data.ipc.enReqState = REQ_STATE_ACK_PENDING;
		return 0;
	}
	case IPC_SET_SHOWIMG_MODE_EVT:
		/* Read the option from the address space of the CGI. */
		bShowThreshold = *((bool*)data.ipc.req.pAddr);
		if (bShowThreshold == TRUE)
		{
			/* Need to capture colored images from now on, this is done in the showRaw state.  */
			STATE_TRAN(me, &me->showThreshold);
		}
		data.ipc.enReqState = REQ_STATE_ACK_PENDING;
		return 0;
	}
	return msg;
}

void MainStateConstruct(MainState *me)
{
	HsmCtor((Hsm *)me, "MainState", (EvtHndlr)MainState_top);
	StateCtor(&me->showRaw, "Show Raw", &((Hsm *)me)->top, (EvtHndlr)MainState_ShowRaw);
	StateCtor(&me->showThreshold, "Show Color", &((Hsm *)me)->top, (EvtHndlr)MainState_ShowThreshold);
}

OscFunction( StateControl)

	OSC_ERR camErr;
	MainState mainState;
	uint8 *pCurRawImg = NULL;

	/* Setup main state machine */
	MainStateConstruct(&mainState);
	HsmOnStart((Hsm *)&mainState);

	OscSimInitialize();

	/* Prologue: initial acquisition setup */
	OscCall( OscCamSetupCapture, OSC_CAM_MULTI_BUFFER);
	OscCall( OscGpioTriggerImage);

	/* Body: infinite acquisition loop */
	while (TRUE)
	{
		/* Wait for captured picture. While a timeout is reported we do service
		 * web interface and GPIO meanwhile. Otherwise we quick. One IPC request
		 * is processed at least. */
		while (TRUE)
		{
			OscCall( HandleIpcRequests, &mainState);

			camErr = OscCamReadPicture(OSC_CAM_MULTI_BUFFER, &pCurRawImg, 0, 4);
			if( camErr == -ETIMEOUT)
			{
				OscCall( HandleIpcRequests, &mainState);
			}
			else
			{
				break;
			}

		}

		/* A valid image is expected. */
		OscAssert_s( camErr == SUCCESS);
		data.pCurRawImg = pCurRawImg;

		/* Process frame by state engine. Sequentially with next capture */
		ThrowEvent(&mainState, FRAMESEQ_EVT);

		/* set new shutter speed */
		if(data.nExposureTimeChanged)
		{
			OscCamSetShutterWidth(data.ipc.state.nExposureTime * 1000);
			data.nExposureTimeChanged = false;
		}

		/* Prepare next capture */
		OscCall( OscCamSetupCapture, OSC_CAM_MULTI_BUFFER);
		OscCall( OscGpioTriggerImage);

		/* Process frame by state engine. Parallel with next capture */
		ThrowEvent(&mainState, FRAMEPAR_EVT);

		/* Advance the simulation step counter. */
		OscSimStep();
	} /* end while ever */

OscFunctionCatch()
OscFunctionEnd()
