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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math-68881.h>

#define M_PI 3.14159265358979323846
#define M_1_PI 0.31830988618379067154

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    double x;
    double y;
    double z; // position, also color (r,g,b)

    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec& norm() { return *this = *this * (1/sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator%(Vec &b){ return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
};

struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()

struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    {                     // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[] = {//Scene: radius, position, emission, color, material
   Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF), //Lite
   Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
   Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
   Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.55,.55,.55),DIFF),//Back
   Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
   Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
   Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
   Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(.4,.4,.3)*.999, SPEC),//Mirr
   Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(.8,.7,.95)*.999, REFR),//Glas
   Sphere(10.5,Vec(23,10.5,98),       Vec(),Vec(.6,1,0.7)*.999, REFR),//Glas
   Sphere(8.,Vec(50,8.,108),       Vec(),Vec(1,0.6,0.7)*.999, REFR),//Glas
   Sphere(6.5, Vec(53,6.5,48),       Vec(),Vec(0.3,.4,.4)*.999, SPEC),//Mirr
 };

const int numSpheres = sizeof(spheres)/sizeof(Sphere);

static inline double clamp(double x)
{
    return x<0 ? 0 : x>1 ? 1 : x;
}

static inline int toInt(double x)
{
    return int(pow(clamp(x),1/2.2)*255+.5);
}

static inline bool intersect(const Ray &r, double &t, int &id)
{
    double n=numSpheres, d, inf=t=1e20;
    for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
    return t<inf;
}
int maximal_ray_depth = 100;

// ca. 650bytes per ray depth, 650KB stack required for max ray depth of 1000

Vec radiance_expl(const Ray &r, int depth, unsigned short *Xi,int E=1){
  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return Vec(); // if miss, return black
  const Sphere &obj = spheres[id];        // the hit object
  Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
  double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
  depth++;

    if (depth > maximal_ray_depth)
        return obj.e*E;
  else
  if (depth>10||!p) { // From depth 10 start Russian roulette
       if (erand48(Xi)<p) f=f*(1/p); else return obj.e*E;
    }
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2.0*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
    Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();

    // Loop over any lights
    Vec e;
    for (int i=0; i<numSpheres; i++){
      const Sphere &s = spheres[i];
      if (s.e.x<=0 && s.e.y<=0 && s.e.z<=0) continue; // skip non-lights

      Vec sw=s.p-x, su=((fabs(sw.x)>.1?Vec(0,1):Vec(1))%sw).norm(), sv=sw%su;
      double cos_a_max = sqrt(1-s.rad*s.rad/(x-s.p).dot(x-s.p));
      double eps1 = erand48(Xi), eps2 = erand48(Xi);
      double cos_a = 1-eps1+eps1*cos_a_max;
      double sin_a = sqrt(1-cos_a*cos_a);
      double phi = 2*M_PI*eps2;
      Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
      l.norm();
      if (intersect(Ray(x,l), t, id) && id==i){  // shadow ray
        double omega = 2*M_PI*(1-cos_a_max);
        e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;  // 1/pi for brdf
      }
    }

    return obj.e*E+e+f.mult(radiance_expl(Ray(x,d),depth,Xi,0));
  } else if (obj.refl == SPEC)              // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance_expl(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
  bool into = n.dot(nl)>0;                // Ray from outside going in?
  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + f.mult(radiance_expl(reflRay,depth,Xi));
  Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
  return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
    radiance_expl(reflRay,depth,Xi)*RP:radiance_expl(Ray(x,tdir),depth,Xi)*TP) :
    radiance_expl(reflRay,depth,Xi)*Re+radiance_expl(Ray(x,tdir),depth,Xi)*Tr);
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t; // distance to intersection
    int id=0; // id of intersected object

    if (!intersect(r, t, id))
        return Vec(); // if miss, return black

    const Sphere &obj = spheres[id];        // the hit object

    Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;

    depth++;

    // Above maximal_ray_depth break recursive loop unconditionally
    if (depth > maximal_ray_depth)
        return obj.e;
    else
    if (depth>5) // From depth of 5 start Russian roulette
    {
        double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl

        if (erand48(Xi)<p)
            f=f*(1/p);
        else return obj.e; //R.R.
    }
    if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
        double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
        Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
    } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
    Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl)>0;                // Ray from outside going in?
    double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj.e + f.mult(radiance(reflRay,depth,Xi));
    Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
        radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
        radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

extern ULONG *workBuffer;
static const int TILE_SIZE = 32;
static const int render_width = TILE_SIZE * 10;
static const int render_height = TILE_SIZE * 8;
static const int samps = 10;
static const int explicit_mode = 1;
int running = TRUE;
struct BitMap *outputBMap = NULL;
struct Window *displayWin = NULL;
struct RastPort *outBMRastPort = NULL;
char tmpbuf[512];
struct timeval start_time;

int GetWidth()
{
    return render_width;
}

int GetHeight()
{
    return render_height;
}

void RedrawTile(int tile_x, int tile_y, int ylines=TILE_SIZE)
{
    struct timeval now;
    WritePixelArray(workBuffer, tile_x * TILE_SIZE, tile_y * TILE_SIZE, render_width * sizeof(ULONG), outBMRastPort,
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
}

void RenderTile(int tile_x, int tile_y)
{
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec cx = Vec(render_width * .5135 / render_height), cy = (cx % cam.d).norm() * .5135, r;

    for (int _y=tile_y * 32; _y < (tile_y + 1) * 32; _y++)
    {
        int y = render_height - _y - 1;
        int render_pos = _y*render_width + tile_x * 32;

        for (unsigned short _x=tile_x * 32, Xi[3]={0,0,(unsigned short)(y*y*y)}; _x < (tile_x + 1) * 32; _x++)   // Loop cols
        {
            int x = _x;
            Vec c = Vec();

            for (unsigned sy=0; sy<2; sy++)     // 2x2 subpixel rows
            {
                for (unsigned sx=0; sx<2; sx++, r=Vec())
                {        // 2x2 subpixel cols
                    for (unsigned s=0; s<samps; s++)
                    {
                        double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/render_width - .5) +
                                cy*( ( (sy+.5 + dy)/2 + y)/render_height - .5) + cam.d;

                        if (explicit_mode)
                            r = r + radiance_expl(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
                        else
                            r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c = c + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
            }

            workBuffer[render_pos++] = ((toInt(c.x) & 0xff)) |
                    ((toInt(c.y) & 0xff) << 8) | ((toInt(c.z) & 0xff) << 16) | 0xff000000;
        }

        RedrawTile(tile_x, tile_y, _y - tile_y * 32 + 1);

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
                        render_width, render_height, 0xC0);
                    EndRefresh(msg->IDCMPWindow, TRUE);
                    break;
            }
            ReplyMsg((struct Message *)msg);
        }
    }
}

void SmallPT(struct Window *myWindow)
{
    struct timeval now;

    displayWin = myWindow;
    running = TRUE;

    outputBMap = AllocBitMap(
                    render_width,
                    render_height,
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

        for (unsigned i=0; i < render_width*render_height; i++)
            workBuffer[i] = 0;

        GetSysTime(&start_time);

        if (explicit_mode)
            spheres[0] = Sphere(5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*20,  Vec(), DIFF);
        else
            spheres[0] = Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF);

        for (unsigned tile_y = 0; tile_y < (render_height / TILE_SIZE); tile_y++)
        {
            for (unsigned tile_x = 0; tile_x < (render_width / TILE_SIZE); tile_x++)
            {
                RenderTile(tile_x, tile_y);
                
                if (running == FALSE)
                    break;
            }
            if (running == FALSE)
                break;
        }
        
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
                        render_width, render_height, 0xC0);
                    EndRefresh(msg->IDCMPWindow, TRUE);
                }

                ReplyMsg((struct Message *)msg);
            }

        } while(running == TRUE);    

        FreeMem(outBMRastPort, sizeof(struct RastPort));
        FreeBitMap(outputBMap);
    }
}