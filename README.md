![alt tag](https://raw.githubusercontent.com/azer89/GPU_Curve_Rendering/master/bezier_large.png)

An implementation of a SIGGRAPH paper:<br/>
Resolution Independent Curve Rendering using Programmable Graphics Hardware (SIGGRAPH 2005)<br/>
by Charles Loop and Jim Blinn<br/>
http://dx.doi.org/10.1145/1073204.1073303<br/>

You can take a look on [CurveRenderer.cpp](https://github.com/azer89/GPU_Curve_Rendering/blob/master/QtTestShader/CurveRenderer.cpp) where the actual algorithm is implemented.

This repo was actually a small part of a bigger project about vectorizing rasterized manga artwork. Below is the video demo:

https://github.com/azer89/GPU_Curve_Rendering/assets/790432/6ef8e14e-4adc-470d-8dd7-3b47f6bafe4e

Dependencies:<br/>
1. Qt 5.x and OpenGL<br/>
2. Visual Studio 2010 <br/>
3. [Eigen Matrix Library](http://eigen.tuxfamily.org)<br/>
