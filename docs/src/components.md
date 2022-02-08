```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using FFTLog
using LaTeXStrings

plot_font = "Computer Modern"
Plots.default(titlefont = (16, plot_font), fontfamily=plot_font,
        linewidth=2, framestyle=:box, fg_legend =:black, label=nothing, grid=false,
        tickfontsize=12, legendfontsize=12, size = (550, 400), labelfontsize = 13,
        dpi = 200)
```


# The FFTLog method
In this page we are going to give the details about the FFTLog method. For a complete reference, we refer to [Hamilton (2000)](https://arxiv.org/abs/astro-ph/9905191) and [Fang et al. (2019)](https://arxiv.org/abs/1911.11947). This method is particularly useful for logarithmically spaced data and is often employed in Cosmology.

## Our objective
We want to perform integrals such a

```math
F(y)=\int_{0}^{\infty} \frac{d x}{x} f(x) j_{\ell}(x y)
```

where ``f(x)`` is the (logarithmically sampled) function we want to transform, ``j_{\ell}`` is the is the order-``\ell`` spherical Bessel function of the first kind.

In the reminder of this documentation, we are gonna show the implementation details.

## The implementation details

The essential idea is to expand ``f(x)`` into a series of power-laws and solve each component integral analytically. We require a logarithmic sampling of ``x`` with linear spacing in ``\ln(x)`` equal to ``\Delta\ln x``, i.e., ``x_q = x_0 \exp(q\Delta\ln x)`` with ``x_0`` being the smallest value in the ``x`` array. The power-law decomposition then means
```math
f\left(x_{q}\right)=\frac{1}{N} \sum_{m=-N / 2}^{N / 2} c_{m} x_{0}^{\nu}\left(\frac{x_{q}}{x_{0}}\right)^{\nu+i \eta_{m}}
```
where ``N`` is the sample size of the input function, ``\eta_m = 2\pi m/(N\Delta\ln x)``, and ``\nu`` is the bias index. The Fourier coefficients
satisfy ``c^∗_m = c_{−m}`` since function ``f(x)`` is real, and are computed by discrete Fourier transforming the “biased” input function
``f(x)/x^\nu`` as

```math
c_{m}=W_{m} \sum_{q=0}^{N-1} \frac{f\left(x_{q}\right)}{x_{q}^{\nu}} e^{-2 \pi i m q / N}
```
where ``W_m`` is a window function which smooths the edges of the ``c_m`` array and takes the form

```math
W(x)= \begin{cases}\frac{x-x_{\min }}{x_{\text {left }}-x_{\min }}-\frac{1}{2 \pi} \sin \left(2 \pi \frac{x-x_{\min }}{x_{\text {left }}-x_{\min }},\right) & x<x_{\text {left }}, \\ 1 & x_{\text {left }}<x<x_{\text {right }}, \\ \frac{x_{\max }-x}{x_{\max }-x_{\text {right }}}-\frac{1}{2 \pi} \sin \left(2 \pi \frac{x_{\max }-x}{x_{\max }-x_{\text {right }}}\right) & x>x_{\text {right }}\end{cases}
```
This is implemented in the following function:

```@docs
FFTLog._c_window
```

```@example tutorial
x = Array(0:1000)
y = FFTLog._c_window(x, 250)
plot(x,y, ylabel = L"W(x)", xlabel = L"x")
```
This filtering is found to reduce the ringing effects.

After the decomposition, the integral reads
```math
F(y)=\frac{1}{N} \sum_{m=-N / 2}^{N / 2} c_{m} x_{0}^{\nu} \int_{0}^{\infty} \frac{d x}{x} \left(\frac{x_{q}}{x_{0}}\right)^{\nu+i \eta_{m}} j_{\ell}(x y)
```

Each term is now analytically solvable. 

```math
\begin{aligned}
F(y) &=\frac{1}{N y^{\nu}} \sum_{m=-N / 2}^{N / 2} c_{m} x_{0}^{-i \eta_{m}} y^{-i \eta_{m}} \int_{0}^{\infty} \frac{d x}{x} x^{\nu+i \eta_{m}} j_{\ell}(x) = \\
&=\frac{\sqrt{\pi}}{4 N y^{\nu}} \sum_{m=-N / 2}^{N / 2} c_{m} x_{0}^{-i \eta_{m}} y^{-i \eta_{m}} g_{\ell}\left(\nu+i \eta_{m}\right)
\end{aligned}
```

where ``g(z)`` is given by

```math
g_{\ell}(z)=2^{z} \frac{\Gamma\left(\frac{\ell+z}{2}\right)}{\Gamma\left(\frac{3+\ell-z}{2}\right)}, \quad-\ell<\Re(z)<2
```

giving the range of bias index ``-\ell<\nu<2``.


Finally, assuming that ``y`` is logarithmically sampled with the same linear spacing ``\Delta \ln y = \Delta x`` in ln y, we can write the last summation as

```math
F\left(y_{p}\right)=\frac{\sqrt{\pi}}{4 y_{p}^{\nu}} \operatorname{IFFT}\left[c_{m}^{*}\left(x_{0} y_{0}\right)^{i \eta_{m}} g_{\ell}\left(\nu-i \eta_{m}\right)\right]
```