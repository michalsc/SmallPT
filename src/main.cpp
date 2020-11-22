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
#include <cybergraphics/cybergraphics.h>
#include <utility/tagitem.h>

#undef SysBase
#define SysBase (*(struct ExecBase **)4)

struct GfxBase *GfxBase = NULL;
struct IntuitionBase *IntuitionBase = NULL;
struct Library *CyberGfxBase = NULL;
struct TimerBase *TimerBase = NULL;
const char * version = VERSION_STRING;
int GetWidth();
int GetHeight();

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
                    TimerBase = (struct TimerBase *)tr->tr_node.io_Device;

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