/* -----------------------------------------
 *              GNUPLOTC
 *
 *  Description: C Interface to gnuplot.
 *  Author: Fernando Muñoz.
 *  Github: ferxinii/gnuplotC
 *  License: MIT
 * --------------------------------------- */

#include "gnuplotc.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#ifdef GNUPLOTC_USE_OMP 
    #include <omp.h>
#endif


// Parallel frame handling
typedef struct frame {
    char *data;
    size_t size;
    struct frame *next;
} s_frame;


typedef struct flist {
    int N;
    s_frame *head;  // Only NOT NULL if in parallel processing mode!
    s_frame *tail;
} s_flist;


void parallel_frame_add(s_flist *flist) {
    s_frame *new_frame = malloc(sizeof(s_frame));
    new_frame->data = NULL;
    new_frame->size = 0;
    new_frame->next = NULL;

    if (flist->head == NULL) {
        flist->head = new_frame;
        flist->tail = new_frame;
        flist->N = 1;  // This should be guaranteed!!
    } else {
        flist->tail->next = new_frame;
        flist->tail = new_frame;
        flist->N++;
    }
}


void free_frame_memory(s_flist *flist) {
    s_frame *current = flist->head;
    while (current) {
        s_frame *next = current->next; // store next before freeing
        free(current->data);
        free(current);
        current = next;
    }
}


// Main interface
typedef struct s_gnuplot {
    FILE *pipe;
    enum gnuplot_type type;
    int dim;
    char *file_name;
    int font_size;
    int size[2];
    int state;
    s_flist frames;
} s_gnuplot;


int GNUPLOTC_FRAMERATE = 24;

int GNUPLOTC_PARAL_VIDEO_THREADS = 0;  


void remove_directory(const char *path) 
{
    DIR *dir = opendir(path);
    if (!dir) return;  // Directory does not exist

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) continue;  // Skip . and ..

        char buf[512];
        snprintf(buf, 512, "%s/%s", path, entry->d_name);

        struct stat statbuf;
        if (stat(buf, &statbuf) == 0) {
            // if (S_ISDIR(statbuf.st_mode)) {
            //     r2 = remove_directory(buf);  // recursive call
            // } 
            if (S_ISREG(statbuf.st_mode)) {  // Test if regular file
                if(unlink(buf) != 0) { // remove file
                    printf("Error: when removing file %s/%s\n", path, entry->d_name);
                }
            } else {
                printf("Error: while removing %s from directory %s. This directory should only contain regular files produced by the code itself.\n", path, entry->d_name);
                exit(1);
            }
        }
    }
    closedir(dir);

    if (rmdir(path) != 0) {  // Directory must be empty before removing it!
        puts("Error: removing directory");
        exit(1);
    }
}


FILE *popen_gnuplot(void)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    if (!pipe) {
        puts("popen_gnuplot: Error.");
        exit(1);
    }
    return pipe;
}


void pipe_terminal_output(s_gnuplot *interface, const char *terminal, char *outfile)
{
    int NMAX = 1024;
    char buff[NMAX];
    char template[] = "set terminal %s enhanced font 'Arial,%d' size %d,%d enhanced\n set output '%s'\n";
    // char template[] = "set terminal %s color font 'Arial,%d' size %d,%d enhanced\n set output '%s'\n";
    int N = -1;

    N = snprintf(buff, NMAX, template, terminal, interface->font_size, 
                 interface->size[0], interface->size[1], outfile);

    if (N < 0 || N >= NMAX) {
        puts("gnuplot_plot: Error with snprintf.");
        exit(1);
    }

    fputs(buff, interface->pipe);
}


void pipe_command(s_gnuplot *interface, char *cmd)
{   // Adds newline \n by default
    if (interface->state != 0) {
        puts("Error! Cannot execute command when plotting is active. Execute all commands before drawing any element.");
        return;
    }
    fputs(cmd, interface->pipe);
    fputs("\n", interface->pipe);
}


void pipe_config_vargs(s_gnuplot *interface, va_list ap)
{   // Pipes configuration variable arguments which may be literal strings or an 
    // array of literal strings (NULL terminated, and marked with a sentinel 
    // GNUPLOTC_ARRAY_MARKER). These can be in any order and as many as needed, but 
    // the last argument must always be NULL.
    char *cmd = va_arg(ap, char *);
    while(cmd != NULL) {
        if (cmd == GNUPLOTC_ARRAY_MARKER) {
            char **cmd_array = va_arg(ap, char **);
            for (int ii=0; cmd_array[ii] != NULL; ii++) {
                pipe_command(interface, cmd_array[ii]);
            }
        } else {
            pipe_command(interface, cmd);
        }
        cmd = va_arg(ap, char *);
    }
}


void gnuplot_config_impl(s_gnuplot *interface, ...)
{
    va_list ap;
    va_start(ap, interface);
    pipe_config_vargs(interface, ap);
    va_end(ap);
}


void pipe_default_commands(s_gnuplot *interface)
{
    pipe_command(interface, "set parametric");
}


void activate_parallel_video_processing(int num_threads)
{
    if (num_threads <= 0) {
        puts("Error: num_threads must be > 0");
        exit(1);
    }
    GNUPLOTC_PARAL_VIDEO_THREADS = num_threads;

#ifdef GNUPLOTC_USE_OMP
    omp_set_num_threads(num_threads);
#endif

    remove_directory(FRAMES_DIR);
    mkdir(FRAMES_DIR, 0755);
}


void setup_video_output(s_gnuplot *interface)
{
    if (GNUPLOTC_PARAL_VIDEO_THREADS == 0) {
        // We pipe each frame from gnuplot directly into ffmpeg as a png
        interface->pipe = popen_gnuplot();

        char buff[256];
        snprintf(buff, 256, "| ffmpeg -y -loglevel error -f png_pipe -r %d -s:v %dx%d -i pipe: -pix_fmt yuv420p -c:v libx264 -crf 18 %s", GNUPLOTC_FRAMERATE, interface->size[0], interface->size[1], interface->file_name);
        pipe_terminal_output(interface, "pngcairo", buff);

    } else if (GNUPLOTC_PARAL_VIDEO_THREADS > 0) {
        // We store each frame as the set of commands to be executed.
        // In the end these templates will be executed in parallel to produce each frame as a png
        // ffmpeg will then make a video from them
        
        parallel_frame_add(&interface->frames);

        interface->pipe = open_memstream(&interface->frames.tail->data, &interface->frames.tail->size);

        char pngfile[256];
        snprintf(pngfile, 256, "%s/frame_%04d.png", FRAMES_DIR, interface->frames.N);
        pipe_terminal_output(interface, "pngcairo", pngfile);
    } else {
        puts("Error: GNUPLOTC_PARAL_VIDEO_THREADS < 0");
        exit(1);
    }
}


s_gnuplot *gnuplot_start_impl(enum gnuplot_type type, char *file_name, int size[2], int font_size, ...) 
{
    s_gnuplot *interface = malloc(sizeof(s_gnuplot));
    interface->pipe = NULL;  // safeguard 
    interface->type = type;
    interface->dim = -1;  // safeguard 
    interface->file_name = strdup(file_name);
    interface->size[0] = size[0];
    interface->size[1] = size[1];
    interface->font_size = font_size;
    interface->state = 0;
    interface->frames.head = NULL;
    interface->frames.N = 1;

    switch (type) {
        case PNG_2D: {
            interface->dim = 2;
            interface->pipe = popen_gnuplot();
            pipe_terminal_output(interface, "pngcairo", file_name);
            break;
        }
        case PDF_2D: {
            interface->dim = 2;
            interface->pipe = popen_gnuplot();
            pipe_terminal_output(interface, "pdfcairo", file_name);
            break;
        }
        case EPS_2D: {
            interface->dim = 2;
            interface->pipe = popen_gnuplot();
            puts("TODO EPS");
            exit(1);
            break;
        }
       case PNG_3D: {
            interface->dim = 3;
            interface->pipe = popen_gnuplot();
            pipe_terminal_output(interface, "pngcairo", file_name);
            break;
        }
        case PDF_3D: {
            interface->dim = 3;
            interface->pipe = popen_gnuplot();
            pipe_terminal_output(interface, "pdfcairo", file_name);
            break;
        }
        case VIDEO_2D: {
            interface->dim = 2;
            setup_video_output(interface);
            break;
        }
        case VIDEO_3D: {
            interface->dim = 3;
            setup_video_output(interface);
            break;
        }
        default: {
            puts("UNIMPLEMENTED TYPE");
            exit(1);
            break;
        }
    }
    
    // Config:
    pipe_default_commands(interface);

    va_list ap;
    va_start(ap, font_size);
    pipe_config_vargs(interface, ap);
    va_end(ap);

    return interface;
}


void next_subplot_impl(s_gnuplot *interface, ...)
{
    interface->state = 0;
    fputs("\n", interface->pipe);
    fflush(interface->pipe);

    va_list ap;
    va_start(ap, interface);
    pipe_config_vargs(interface, ap);
    va_end(ap);
}


void next_frame_impl(s_gnuplot *interface, ...)
{
    interface->state = 0;
    fputs("\n", interface->pipe);
    fflush(interface->pipe);

    if (GNUPLOTC_PARAL_VIDEO_THREADS > 0) {
        // Open new file to write template of frame to:
        fclose(interface->pipe);
        
        parallel_frame_add(&interface->frames);
        interface->pipe = open_memstream(&interface->frames.tail->data, &interface->frames.tail->size);

        // Setup config
        char buff[256];
        snprintf(buff, 256, "%s/frame_%04d.png", FRAMES_DIR, interface->frames.N);
        pipe_terminal_output(interface, "pngcairo", buff);
        pipe_default_commands(interface);
    } else {
        interface->frames.N++;
    }
    
    va_list ap;
    va_start(ap, interface);
    pipe_config_vargs(interface, ap);
    va_end(ap);
}



void process_video_parallel(s_gnuplot *interface)
{   
    // Create array of pointers to frame for easy reference
    s_frame **farray = malloc(sizeof(s_frame*) * interface->frames.N);
    s_frame *current = interface->frames.head;
    for (int ii=0; ii<interface->frames.N; ii++) {
        farray[ii] = current;
        current = current->next;
    }

#ifndef GNUPLOTC_USE_OMP 
    int active_procs = 0;
    for (int ii=0; ii<interface->frames.N; ii++) {
        // Limit number of concurrent processes
        while (active_procs >= GNUPLOTC_PARAL_VIDEO_THREADS) {
            wait(NULL);  // Wait for any child
            active_procs--;
        }

        pid_t pid = fork();
        if (pid < 0) { perror("fork"); exit(1); } 
        else if (pid == 0) {  // Child process
            FILE *gp = popen("gnuplot", "w");
            if (!gp) { perror("popen gnuplot"); exit(1); }
            fwrite(farray[ii]->data, 1, farray[ii]->size, gp);
            pclose(gp);  // Wait for gnuplot to finish
            exit(0);      // Child exits
        } 
        else {  // Parent process: track active children
            active_procs++;
        }
    }
    while (active_procs > 0) {  // Wait for remaining children
        wait(NULL);
        active_procs--;
    }
#else
    if (GNUPLOTC_PARAL_VIDEO_THREADS == 1) {
        for (int ii=0; ii<interface->frames.N; ii++) {
            FILE *gp = popen("gnuplot","w");
            fwrite(farray[ii]->data, 1, farray[ii]->size, gp);
            pclose(gp);
        }
    } else {
        #pragma omp parallel for schedule(dynamic)
        for (int ii=0; ii<interface->frames.N; ii++) {
            s_frame *frame = farray[ii];
            FILE *gnuplot_instance = popen("gnuplot", "w");
            fwrite(frame->data, 1, frame->size, gnuplot_instance);
            pclose(gnuplot_instance);
        }
    }
#endif

    free_frame_memory(&interface->frames);
    free(farray);
    
    char buff[1024];
    snprintf(buff, 1024, "ffmpeg -y -loglevel error -framerate %d -s:v %dx%d -i %s/frame_%%04d.png -pix_fmt yuv420p -c:v libx264 -crf 18 %s", GNUPLOTC_FRAMERATE, interface->size[0], interface->size[1], FRAMES_DIR, interface->file_name);
    system(buff);
}



void gnuplot_end(s_gnuplot *interface)
{   
    fputs("\n", interface->pipe);
    if (interface->type != VIDEO_3D && interface->type != VIDEO_2D) {
        pclose(interface->pipe);
    } else {
        if (GNUPLOTC_PARAL_VIDEO_THREADS == 0) {
            pclose(interface->pipe);
        } else if (GNUPLOTC_PARAL_VIDEO_THREADS > 0) {
            fclose(interface->pipe);
            process_video_parallel(interface);
            remove_directory(FRAMES_DIR);
        } else {
            puts("Error: GNUPLOTC_PARAL_VIDEO_THREADS < 0.");
        }
    }
    free(interface->file_name);
}


void video_to_gif(const char *file_video, const char *file_gif, int size[2], int framerate)
{
      char command_mp4_to_gif[1024];
      snprintf(command_mp4_to_gif, 1024, "ffmpeg -loglevel error -i %s -vf \"fps=%d,scale=%d:%d:flags=lanczos\" -gifflags +transdiff -y %s 2>&1", file_video, framerate, size[0], size[1], file_gif);
      system(command_mp4_to_gif);
}


int guard_is2d(s_gnuplot *interface)
{
    if (interface->dim == 2) return 1;
    puts("Error: Cannot use 2d method!");
    return 0;
}


int guard_is3d(s_gnuplot *interface)
{
    if (interface->dim == 3) return 1;
    puts("Error: Cannot use 3d method!");
    return 0;
}


void guard_active_plotting(s_gnuplot *interface)
{
    if (interface->state == 0) {
        if (interface->dim == 2) {
            fputs( "plot ", interface->pipe);
        } else if (interface->dim == 3) {
            fputs( "splot ", interface->pipe);
        } else {
            puts("Error. dim not valid."); 
            exit(1);
        }
        interface->state = 1;
    } else if (interface-> state == 1) {
        fputs(" , ", interface->pipe);
    } else {
        puts("Error. State is not valid.");
    }
}


void add_datablock_from_file(s_gnuplot *interface, const char *filename, int header_lines_skip, const char *datablock_name)
{
    if (interface->state != 0) {
        puts("Cannot add datablock when plotting is active. Do so before drawing any element."); 
        return;
    }

    FILE *file = fopen(filename, "r");
    if (!file) {puts("Error opening file for data block in plot"); return;}

    fprintf(interface->pipe, "$%s << EOD\n", datablock_name);
    char buffer[1024];

    for (int ii=0; ii<header_lines_skip; ii++) {
        fgets(buffer, sizeof(buffer), file);
    }
    while (fgets(buffer, sizeof(buffer), file)) {
        fputs(buffer, interface->pipe);
    }
    fputs("EOD\n", interface->pipe);
    fclose(file);
}


int config_specifies_title(const char *config)
{
    return strstr(config, "title") != NULL || strstr(config, "notitle") != NULL;
}


void pipe_element_type(s_gnuplot *interface, enum element_type type)
{
    fprintf(interface->pipe, " w ");
    switch (type) {
        case LINES: {
            fprintf(interface->pipe, "lines");
            break;
        } case POINTS: {
            fprintf(interface->pipe, "points");
            break;
        } case LINESPOINTS: {
            fprintf(interface->pipe, "linespoints");
            break;
        } case POLYGONS: {
            fprintf(interface->pipe, "polygons");
            break;
        } case PM3D: {
            fprintf(interface->pipe, "pm3d");
            break;
        } default: {
            puts("Error: Element type not implemented!");
            exit(1);
        }
    }
}

void pipe_element_config(s_gnuplot *interface, const char *config)
{
    if (!config) {
        fprintf(interface->pipe, " notitle");
        return;
    }

    if (config_specifies_title(config)) fprintf(interface->pipe, " %s", config);
    else fprintf(interface->pipe, " %s notitle", config);
}


void draw_datablock(s_gnuplot *interface, const char *datablock_name, enum element_type type, const char *config)
{
    guard_active_plotting(interface);
    fprintf(interface->pipe, "$%s", datablock_name);
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_file(s_gnuplot *interface, const char *file_name, int header_lines_skip, enum element_type type, const char *config)
{
    guard_active_plotting(interface);
    fprintf(interface->pipe, "'%s' skip %d", file_name, header_lines_skip);
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_command(s_gnuplot *interface, const char *command)
{
    guard_active_plotting(interface);
    fprintf(interface->pipe, "%s ", command);
}


void draw_point_2d(s_gnuplot *interface, double x, double y, const char *config)
{
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'%f %f\'\" w p", x, y);
    pipe_element_config(interface, config);
} 


void draw_2d(s_gnuplot *interface, double *x, double *y, int N, enum element_type type, const char *config)
{
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        fprintf(interface->pipe, "%f %f", x[ii], y[ii]);
        if (ii < N-1) {
            fprintf(interface->pipe, "\\n");
        }
    }
    fprintf(interface->pipe, "\'\"");
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_array_2d(s_gnuplot *interface, double **coords, int N, enum element_type type, const char *config)
{ 
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        fprintf(interface->pipe, "%f %f", coords[ii][0], coords[ii][1]);
        if (ii < N-1) {
            fprintf(interface->pipe, "\\n");
        }
    }
    fprintf(interface->pipe, "\'\"");
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_errorbars_2d(s_gnuplot *interface, double *x, double *mean, double *err, int N, const char *config)
{
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        fprintf(interface->pipe, "%f %f %f %f", x[ii], mean[ii], mean[ii] + err[ii], mean[ii] - err[ii]);
        if (ii < N-1) {
            fprintf(interface->pipe, "\\n");
        }
    }
    fprintf(interface->pipe, "\'\" w errorbars");
    pipe_element_config(interface, config);
}


void draw_arrows_from_array_2d(s_gnuplot *interface, double **coords, int N, int spacing, int offset, const char *config)
{   
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        if ((offset +ii)%spacing == 0 && ii+1 < N) {
            fprintf(interface->pipe, "%f %f %f %f", coords[ii][0], coords[ii][1], coords[ii+1][0]-coords[ii][0], coords[ii+1][1]-coords[ii][1]);
        }
        if (ii < N-1) {
            fprintf(interface->pipe, "\\n");
        }
    }
    fprintf(interface->pipe, "\'\" w vectors");
    pipe_element_config(interface, config);
}


void draw_segment_2d(s_gnuplot *interface, double x0, double y0, double xf, double yf, const char *config)
{
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'%f %f %f %f\' \" w vectors nohead", x0, y0, xf-x0, yf-y0);
    pipe_element_config(interface, config);
}


void draw_function_2d(s_gnuplot *interface, double x0, double xf, int N, double (*fun)(double), enum element_type type, const char *config)
{
    if (!guard_is2d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N; ii++) {
        double x = x0 + ii * (xf - x0) / (N - 1);
        fprintf(interface->pipe, "%f %f", x, fun(x));
        if (ii < N-1) {
            fprintf(interface->pipe, "\\n");
        }
    }
    fprintf(interface->pipe, "\'\"");
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_sphere_3d(s_gnuplot *interface, double x, double y, double z, double r, enum element_type type, const char *config)
{   // type must be: pm3d, lines, points, ... (polygons? Should produce nothing if haven't "set pm3d", else, same as pm3d)
    if (!guard_is3d(interface)) return;
    guard_active_plotting(interface);
    fprintf(interface->pipe, "[u=-pi/2:pi/2] [v=0:2*pi] %f+%f*cos(u)*cos(v),%f+%f*cos(u)*sin(v),%f+%f*sin(u)", 
            x, r, y, r, z, r);
    pipe_element_type(interface, type);
    pipe_element_config(interface, config);
}


void draw_solid_triangle_3d(s_gnuplot *interface, double v0[3], double v1[3], double v2[3], const char *config)
{
    if (!guard_is3d(interface)) return;
    guard_active_plotting(interface);
    fprintf(interface->pipe, "\"<echo \'");
    fprintf(interface->pipe, "%f %f %f\\n", v0[0], v0[1], v0[2]);
    fprintf(interface->pipe, "%f %f %f\\n", v1[0], v1[1], v1[2]);
    fprintf(interface->pipe, "%f %f %f\\n\'\"", v2[0], v2[1], v2[2]);
    fprintf(interface->pipe, "w polygons");
    pipe_element_config(interface, config);
}


void draw_polytope_3d(s_gnuplot *interface, double **coords, int N, const char *config)
{
    if (!guard_is3d(interface)) return;
    guard_active_plotting(interface);

    fprintf(interface->pipe, "\"<echo \'");
    for (int ii=0; ii<N-1; ii++) {
        fprintf(interface->pipe, "%f %f %f\\n", coords[ii][0], coords[ii][1], coords[ii][2]);
    }
    fprintf(interface->pipe, "%f %f %f\'\"", coords[N-1][0], coords[N-1][1], coords[N-1][2]);
    fprintf(interface->pipe, "w polygons");
    pipe_element_config(interface, config);
}


