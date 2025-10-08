# gnuplotC

Simple and easy C interface to *gnuplot*. Allows to draw figures and create video animations (efficient by using parallelism).

Currently, **gnuplotC** supports the following modes:


|    | **png** | **pdf** | **eps** | **video (mp4)** |
| ---: | :---: | :---: | :---: | :---: |
| **2D** | PNG_2D | PDF_2D | EPS_2D | VIDEO_2D |
| **3D** | PNG_3D | PDF_3D |  | VIDEO_3D |




---
### Creating a plot

1. **Start *gnuplot* interface.** Select mode, file name, figure size, ...
2. *(Optional)* **Set *gnuplot* configurations.** Before adding elements!
3. **Add elements.**
4. **End.**

```C
int fontsize = 18, figsize[2] = {1080, 1080};

t_gnuplot *ifc = gnuplot_start(PNG_2D, "segment.png", figsize, fontsize);
draw_segment_2d(ifc, 1, 0, 0.4, 0.6, "lw 2"); 
gnuplot_end(ifc);
```



---
### Creating a video
Videos are built using *ffmpeg*, which must be installed and available. Treat each frame as its own independent figure. By default, each frame is piped directly to *ffmpeg* without the need to store it.

1. **Start *gnuplot* interface.** Select framerate, ...
2. *(Optional)*  **Set configurations.** 
3. **Add elements.**
4. **Start a new frame.**
5. **Repeat from step 2 as needed.**
6. **End.**

Optionally, the default framerate (24 fps) can be changed by modifying the global variable *GNUPLOTC_FRAMERATE* **before** starting the interface.

```C
GNUPLOTC_FRAMERATE = 12;

char *video_config[] = {"unset border", "unset tics", "set isosamples 100", 
                        "set view equal xyz", NULL};

ifc = gnuplot_start(VIDEO_3D, "sphere.mp4", figsize, fontsize, 
                    GNUPLOTC_ARRAY(video_config), "set view 70, 0");
draw_sphere_3d(ifc, 0.5, 0.5, 0.5, 0.2, LINES, NULL);
for (int ii=1; ii<GNUPLOTC_FRAMERATE*3+1; ii++) {
    char buff[256]; 
    snprintf(buff, 256, "set view 70, %f", fmod(ii*360.0/(GNUPLOTC_FRAMERATE*3), 360.0));
    next_frame(ifc, GNUPLOTC_ARRAY(video_config), buff);
    draw_sphere_3d(ifc, 0.5, 0.5, 0.5, 0.2, LINES, NULL);
}
gnuplot_end(ifc);
```



### Faster parallel video processing
Parallel processing of the frames is possible, with each thread running *gnuplot* in parallel resulting in much faster processing times. The downside is that all the frames must be saved to disc beforehand (in a temporal directory), which can take up some space.

To activate this mode, simply call *activate_parallel_video_processing(...)* **before** starting the interface.



---
### Configuring the interface

The gnuplot configuration can be specified:
- When calling *gnuplot_start(...)*, in the form of additional optional parameters.
- Calling *gnuplot_config(...)*, **before** drawing any element.
- In the case of a video, also when calling *new_frame(...)*.

Once an element in a plot (or frame) is drawn, it is not possible to change the gnuplot configuration.

```C
t_gnuplot *ifc = gnuplot_start(PNG_2D, "circle.png", figsize, fontsize, "set title 'CIRCLE'", "set xlabel 'x'");
gnuplot_config(ifc, "set xrange [0:1]", "set yrange [0:1]");
[...]
```


Also, we may define a **NULL terminated** array of configurations, which we can pass at any point as an additional configuration using *GNUPLOTC_ARRAY()*:
```C
char *cmd_array[] = {"set xrange [0:1]", "set yrange [0:1]", NULL};
t_gnuplot *ifc = gnuplot_start(PNG_2D, "test.png", figsize, fontsize, GNUPLOTC_ARRAY(cmd_array));
[...]
```



--- 
### Dependencies and compilation

The user must have available *gnuplot* and *ffmpeg* (for video creation) as command line tools. Also, 

#### Optional OMP mode
By default, parallel processing of frames is implemented using *fork()* to spawn child processes. However, those users who prefer an *OMP* approach, may activate it one of two ways:
- Defining *GNUPLOTC_USE_OMP* before including *gnuplotc.h*.
    ```C
    #define GNUPLOTC_USE_OMP
    #include "gnuplotc.h"
    ```
- Or by passing *-DGNUPLOTC_USE_OMP* as an argument to the compiler.
    ```
    gcc-15 -DGNUPLOTC_USE_OMP example.c plot.c -o example_omp
    ```
In this case, the compiler should support OMP (in my Macbook the native gcc does not and I have to use a Homebrew version).


---
### Examples

**Basic examples (see /example/example.c)**

<p align="center">
<img src="./example/10.png" alt="Example 1" width="400" height="auto" />
<img src="./example/16.gif" alt="Example 2" width="400" height="auto">
</p>

**Other examples I made using the library**
<p align="center">
<img src="./example/ex_1.png" alt="Example 3" width="400" height="auto" />
<img src="./example/ex_2.png" alt="Example 4" width="400" height="auto">
</p>
<p align="center">
<img src="./example/ex_3.png" alt="Example 5" width="400" height="auto" />
<img src="./example/ex_3d.png" alt="Example 6" width="400" height="auto">
</p>
