# autocode-SCvx
We address the embedded code generation for an optimal control algorithm, [SCvx](https://arxiv.org/abs/1804.06539) , which is particularly suitable for solving trajectory planning problems with collision avoidance constraints. Existing uses of SCvx on drones or embedded platform are currently custom coded.
On the other hand, recent toolboxes such as [SCPToolbox](https://github.com/UW-ACL/SCPToolbox.jl) provide a simpler access to these trajectory planning algorithms, based on the resolution of a sequence of convex sub-problems.

We define here a framework, in Python, enabling the automatic code generation for SCvx , in C, based on [CVXPYgen](https://github.com/cvxgrp/cvxpygen) and the ECOS solver. The framework is able to address problems involving non- convex constraints such as obstacle avoidance.
This is a first step towards a more streamlined process to autocode trajectory planning algorithms and convex optimization solvers.

## Use case
In the main file, you can find the parameters describing a Dubin's car trajectory optimization problem, with two obstacles.

### Environnement
Using Anaconda, you can just import [scvx_env.yaml](https://github.com/Dberrah/autocode-SCvx/blob/main/scvx_env.yaml) as a new environnement.
To compile and run the C code, ``CMake 3.5`` or newer is required.

### Example
Run the following command to solve the problem and/or automatically generate the custom solver in C.
```bash
python main.py
```

On Unix platforms, run the following commands in your terminal to compile and run the custom C solver:

```bash
cd SCvx/c/build
cmake ..
cmake --build . --target SCvx
./SCvx
```
