# RieszDSP.jl

RieszDSP is a Julia implementation of the forward and inverse transform of the higher-order Riesz-wavelet frame from (Unser et. al. 2011). This transform takes a multi-dimensional array as input, constructs a higher-order Riesz transform, followed by a wavelet transform using a isotropic template filter. So far, only the Simoncelli filter from Unser's article is implemented. This is a perfect reconstruction transform, which means the round-trip discrepancy (i.e., difference between the input array and the array from applying the of the Riesz-wavelet transform and then applying the inverse Riesz-wavelet transform) should be zero, up to numerical precision.

A portion of the code here is based on the [Generalized Riesz-Wavelet Toolbox for Matlab](https://bigwww.epfl.ch/demo/steerable-wavelets/) authored by Nicolas Chenouard, Dimitri Van De Ville and Michael Unser. Their toolbox as well as RieszDSP use a frequency-domain implementation of filtering operations.

The implementation in RieszDSP does not use critical sampling, so each transform response is an array of the same size as the input array. The original implementation from the MATLAB toolbox uses critical sampling, so the transform responses are not arrays of the same size as the input array. See [Demo: image analysis](@ref) for an example.

# Specification of the higher-order Riesz transform

Let the dimension of the input array be ``D``. Given a positive integer ``L``, the ``L``-th order Riesz transform is a collection of iterated Riesz transforms, each of which is specified by a multi-index ``a \in \mathbb{A}\left(D,L\right)``. The index set is defined as:
```math
\begin{align*}
\left[0:L\right] & :=\left\{ 0,1,\cdots,L\right\} ,\\
\left[0:L\right]^{D} & :=\underbrace{\left[0:L\right]\times\left[0:L\right]\times\cdots\times\left[0:L\right]}_{\text{D times}},\\
\mathbb{A}\left(D,L\right) & :=\left\{ a\in\left[0:L\right]^{D}\,\mid\,\left|a\right|=L\right\} \subset\left[0:L\right]^{D}.
\end{align*}
```
The size of the index set is 
```math
\begin{align*}
\left|\mathbb{A}\left(D,L\right)\right| & =\begin{pmatrix}L+D-1\\
D-1
\end{pmatrix}.
\end{align*}
```
See Unser 2011 for further details. The set ``\mathbb{A}\left(D,L\right)`` is denoted by `a_array` in our image analysis demo, and is one of the returned variables from [`rieszwaveletanalysis`](@ref)

We refer to the multi-index ``a`` as the *Riesz order state* throughout this documentation. The contents of the state ``a`` has ``D`` elements, each of which is an integer that takes value from ``0`` to ``L``. These constraints allow us to use the state ``a`` to specify how many times the Riesz transform is applied to the input in the first dimension axis, in the second dimension axis, ..., to the ``D``-th dimension axis. An interesting analogy is that if we replace Riesze transform with the partial derivative operator that acts along one of the ``D`` coordinate axes, then we can use ``a`` to specify a ``L``-th order partial derivative operator. Essentially, the set of states ``\mathbb{A}\left(D,L\right)`` is the index set for a symmetric tensor. The Riesz transform (and the partial derivative operator for an infinitely differentiable function) along an axis is commutative, that is why we end up with a symmetric tensor of iterated Riesz transforms.

In the demo code, you'll see that the Riesz transform is like a smoothed gradient operator. Given a Riesz order state ``a``, you'll see the edges along the axis where the transform was repeated (according to the contents of ``a``) has a higher magnitude.

# Specification of the wavelet analysis/synthesis
The nomenclature used in this documentation is *wavelet analysis* for the decomposition of an input multi-dimensional array to several multi-dimensional arrays, each called a *wavelet subband response* and corresponds to a different *scale*.

Given an input `y`, a heuristic formula for choosing the number of scales of the wavelet analysis (i.e., the number of subband responses we'd get) is:

```math
\begin{align*}
N_{\text{scales}}= & \left\lfloor \log_{2}\left(\max\left\{ N_{1},N_{2},\cdots,N_{D}\right\} \right)\right\rfloor 
\end{align*}
```
 where ``N_{d}`` is the size of the grid in the ``d``-th dimension.
For example, a 256 by 512 grayscale image has ``N_{1}=256``, and ``N_{2}=512``.

# References
1. M. Unser, N. Chenouard, D. Van de Ville, Steerable Pyramids and Tight Wavelet Frames in L2(Rd), IEEE Transactions on Image Processing, 2011. [DOI: 10.1109/TIP.2011.2138147](https://doi.org/10.1109/TIP.2011.2138147).

# Acknowledgement
I thank Nicolas Chenouard and Dimitri Van De Ville for answering questions about their toolbox and the Riesz-wavelet transform.