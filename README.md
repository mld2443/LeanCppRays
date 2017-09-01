# LeanRay
After discussing this, it was decided that I would rewrite it in C++, and that it needs to be lean, and that I do not need to spend a lot of time documenting the code and explaining how everything works in comments and docs. As part of it being lean, I chose to leave everything in one admittedly long file, (an honored C/C++ tradition).

![1920x1080 3000](https://user-images.githubusercontent.com/5340992/29981586-c2fce13c-8f13-11e7-8db2-ead4431d95e7.png)
## Running
This program requires [libpng](http://www.libpng.org/pub/png/libpng.html) to compile. To run the ray tracer, just compile and link the library, then execute it.
If you wish to alter the scene, the scene description is at the bottom of the file in the main function.

### How does my ray tracer work?
This ray tracer is almost a complete port of my previous [Java ray tracer](https://github.com/mld2443/SingleThreadedRayTracer). See that project for more information on how my specific ray tracer works. For more general information on how ray tracing works, [Wikipedia has a fairly user-friendly description](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)).

## Concerns
* Since I wrote this with as little documentation as possible, I was certain to name everything very clearly.
* All parameters are designated at the bottom of the file in the main function. The choice not to have file input or command-line parameters was out of concern for bloating the already long file.
* I initially had issues with memory leakage, using up several gigabytes of space and rising while tracing. This was solved quickly. The tracing is stable, and I made sure that even the program completion is free of leaks.
