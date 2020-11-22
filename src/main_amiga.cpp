#define L 
#include <exec/types.h>
#include <exec/execbase.h>
#include <intuition/intuition.h>
#include <intuition/screens.h>
#include <workbench/startup.h>
#include <graphics/gfxbase.h>
#include <graphics/gfx.h>
#include <dos/dos.h>
#include <dos/dosextens.h>
#include <exec/io.h>
#include <devices/timer.h>
#include <proto/exec.h>
#include <proto/graphics.h>
#include <proto/intuition.h>
#include <proto/cybergraphics.h>
#include <proto/timer.h>
#include <cybergraphics/cybergraphics.h>
#include <utility/tagitem.h>
#include <stdio.h>
#include "smallpt.h"

#undef SysBase
#define SysBase (*(struct ExecBase **)4)

struct GfxBase *GfxBase = NULL;
struct IntuitionBase *IntuitionBase = NULL;
struct Library *CyberGfxBase = NULL;
struct Device *TimerBase = NULL;
const char * version = VERSION_STRING;

struct Window * createMainWindow(int req_width, int req_height)
{
    struct Window *displayWin = NULL;

    displayWin = OpenWindowTags(0,
                                WA_InnerWidth, req_width,
                                WA_InnerHeight, req_height,
                                WA_Title, (ULONG) "SmallPT renderer",
                                WA_SimpleRefresh, TRUE,
                                WA_CloseGadget, TRUE,
                                WA_DepthGadget, TRUE,
                                WA_DragBar, TRUE,
                                WA_SizeGadget, FALSE,
                                WA_SizeBBottom, FALSE,
                                WA_SizeBRight, FALSE,
                                WA_IDCMP, IDCMP_CLOSEWINDOW | IDCMP_REFRESHWINDOW,
                                TAG_DONE);

    return displayWin;
}

ULONG *workBuffer;

void SmallPT(struct Window *myWindow);

int main(void)
{
    IntuitionBase = (struct IntuitionBase *)OpenLibrary((CONST_STRPTR)"intuition.library", 0);
    GfxBase = (struct GfxBase *)OpenLibrary((CONST_STRPTR) "graphics.library", 0);
    CyberGfxBase = OpenLibrary((CONST_STRPTR)"cybergraphics.library", 0);
    struct Window *myWindow = NULL;
    struct MsgPort *timerPort = CreateMsgPort();
    struct timerequest *tr = (struct timerequest *)CreateIORequest(timerPort, sizeof(struct timerequest));

    if (IntuitionBase != NULL && GfxBase != NULL && CyberGfxBase != NULL)
    {
        if (timerPort)
        {
            FreeSignal(timerPort->mp_SigBit);
            timerPort->mp_SigBit = -1;
            timerPort->mp_Flags = PA_IGNORE;

            if (tr)
            {
                if (!OpenDevice("timer.device", UNIT_VBLANK, (struct IORequest *)tr, 0))
                {
                    TimerBase = (struct Device *)tr->tr_node.io_Device;

                    myWindow = createMainWindow(GetWidth(), GetHeight());

                    if (myWindow != NULL)
                    {
                        workBuffer = (ULONG *)AllocMem(GetWidth() * GetHeight() * sizeof(ULONG), MEMF_ANY|MEMF_CLEAR);
                        
                        if (workBuffer != NULL)
                        {
                            SmallPT(myWindow);

                            FreeMem(workBuffer, GetWidth() * GetHeight() * sizeof(ULONG));
                        }

                        CloseWindow(myWindow);
                    }

                    CloseDevice((struct IORequest *)tr);
                }

                DeleteIORequest((struct IORequest *)tr);
            } 

            DeleteMsgPort(timerPort);
        }
    }

    if (IntuitionBase != NULL)
        CloseLibrary((struct Library *)IntuitionBase);

    if (GfxBase != NULL)
        CloseLibrary((struct Library *)GfxBase);

    if (CyberGfxBase != NULL)
        CloseLibrary((struct Library *)CyberGfxBase);

    return 0;
}

struct BitMap *outputBMap = NULL;
struct Window *displayWin = NULL;
struct RastPort *outBMRastPort = NULL;
char tmpbuf[512];
struct timeval start_time;

void RedrawTile(int tile_x, int tile_y, int ylines)
{
    struct timeval now;
    WritePixelArray(workBuffer, tile_x * TILE_SIZE, tile_y * TILE_SIZE, GetWidth() * sizeof(ULONG), outBMRastPort,
                                tile_x * TILE_SIZE, tile_y * TILE_SIZE, TILE_SIZE, ylines, RECTFMT_ARGB);

    BltBitMapRastPort (outputBMap, tile_x * TILE_SIZE, tile_y * TILE_SIZE, displayWin->RPort,
                                displayWin->BorderLeft + tile_x * TILE_SIZE, displayWin->BorderTop + tile_y * TILE_SIZE,
                                            TILE_SIZE, ylines, 0xC0);
    
    GetSysTime(&now);
    SubTime(&now, &start_time);
    snprintf(tmpbuf, 500, "SmallPT renderer: %ld:%02ld:%02ld (working...)", 
                                    now.tv_secs / 3600,
                                    (now.tv_secs / 60) % 60,
                                    now.tv_secs % 60);
    SetWindowTitles(displayWin, tmpbuf, NULL);

    struct IntuiMessage *msg;
    while ((msg = (struct IntuiMessage *)GetMsg(displayWin->UserPort)))
    {
        switch(msg->Class)
        {
            case IDCMP_CLOSEWINDOW:
                running = FALSE;
                return;

            case IDCMP_REFRESHWINDOW:
                BeginRefresh(msg->IDCMPWindow);
                BltBitMapRastPort (outputBMap, 0, 0,
                    msg->IDCMPWindow->RPort, msg->IDCMPWindow->BorderLeft, msg->IDCMPWindow->BorderTop,
                    GetWidth(), GetHeight(), 0xC0);
                EndRefresh(msg->IDCMPWindow, TRUE);
                break;
        }
        ReplyMsg((struct Message *)msg);
    }
}

void SmallPT(struct Window *myWindow)
{
    struct timeval now;

    displayWin = myWindow;
    running = TRUE;

    outputBMap = AllocBitMap(
                    GetWidth(),
                    GetHeight(),
                    GetBitMapAttr(myWindow->WScreen->RastPort.BitMap, BMA_DEPTH),
                    BMF_DISPLAYABLE, myWindow->WScreen->RastPort.BitMap);

    if (outputBMap != NULL)
    {
        outBMRastPort = (struct RastPort *)AllocMem(sizeof(struct RastPort), MEMF_ANY|MEMF_CLEAR);
        if (!outBMRastPort)
        {
            FreeBitMap(outputBMap);
            return;
        }

        InitRastPort(outBMRastPort);
        outBMRastPort->BitMap = outputBMap;

        for (unsigned i=0; i < GetWidth() * GetHeight(); i++)
            workBuffer[i] = 0;

        GetSysTime(&start_time);

        SmallPTLoop();
        
        GetSysTime(&now);
        SubTime(&now, &start_time);
        snprintf(tmpbuf, 500, "SmallPT renderer: Total time %ld:%02ld:%02ld", 
                                        now.tv_secs / 3600,
                                        (now.tv_secs / 60) % 60,
                                        now.tv_secs % 60);
        SetWindowTitles(displayWin, tmpbuf, NULL);
        
        running = TRUE;

        do {
            struct IntuiMessage *msg;

            WaitPort(myWindow->UserPort);
            while ((msg = (struct IntuiMessage *)GetMsg(myWindow->UserPort)) != 0)
            {
                if (msg->Class == IDCMP_CLOSEWINDOW)
                {
                    running = FALSE;
                }
                else if (msg->Class == IDCMP_REFRESHWINDOW)
                {
                    BeginRefresh(msg->IDCMPWindow);
                    BltBitMapRastPort (outputBMap, 0, 0,
                        msg->IDCMPWindow->RPort, msg->IDCMPWindow->BorderLeft, msg->IDCMPWindow->BorderTop,
                        GetWidth(), GetHeight(), 0xC0);
                    EndRefresh(msg->IDCMPWindow, TRUE);
                }

                ReplyMsg((struct Message *)msg);
            }

        } while(running == TRUE);    

        FreeMem(outBMRastPort, sizeof(struct RastPort));
        FreeBitMap(outputBMap);
    }
}
