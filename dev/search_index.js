var documenterSearchIndex = {"docs":
[{"location":"components/","page":"Components","title":"Components","text":"using Plots; gr()\nPlots.reset_defaults()\nusing FFTLog\nusing LaTeXStrings\n\nplot_font = \"Computer Modern\"\nPlots.default(titlefont = (16, plot_font), fontfamily=plot_font,\n        linewidth=2, framestyle=:box, fg_legend =:black, label=nothing, grid=false,\n        tickfontsize=12, legendfontsize=12, size = (550, 400), labelfontsize = 13,\n        dpi = 200)","category":"page"},{"location":"components/#The-FFTLog-method","page":"Components","title":"The FFTLog method","text":"","category":"section"},{"location":"components/","page":"Components","title":"Components","text":"In this page we are going to give the details about the FFTLog method. For a complete reference, we refer to Hamilton (2000) and Fang et al. (2019). This method is particularly useful for logarithmically spaced data and is often employed in Cosmology.","category":"page"},{"location":"components/#Our-objective","page":"Components","title":"Our objective","text":"","category":"section"},{"location":"components/","page":"Components","title":"Components","text":"We want to perform integrals such a","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"F(y)=int_0^infty fracd xx f(x) j_ell(x y)","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"where f(x) is the (logarithmically sampled) function we want to transform, j_ell is the is the order-ell spherical Bessel function of the first kind.","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"In the reminder of this documentation, we are gonna show the implementation details.","category":"page"},{"location":"components/#The-implementation-details","page":"Components","title":"The implementation details","text":"","category":"section"},{"location":"components/","page":"Components","title":"Components","text":"The essential idea is to expand f(x) into a series of power-laws and solve each component integral analytically. We require a logarithmic sampling of x with linear spacing in ln(x) equal to Deltaln x, i.e., x_q = x_0 exp(qDeltaln x) with x_0 being the smallest value in the x array. The power-law decomposition then means","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"fleft(x_qright)=frac1N sum_m=-N  2^N  2 c_m x_0^nuleft(fracx_qx_0right)^nu+i eta_m","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"where N is the sample size of the input function, eta_m = 2pi m(NDeltaln x), and nu is the bias index. The Fourier coefficients satisfy c^_m = c_m since function f(x) is real, and are computed by discrete Fourier transforming the “biased” input function f(x)x^nu as","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"c_m=W_m sum_q=0^N-1 fracfleft(x_qright)x_q^nu e^-2 pi i m q  N","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"where W_m is a window function which smooths the edges of the c_m array and takes the form","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"W(x)= begincasesfracx-x_min x_text left -x_min -frac12 pi sin left(2 pi fracx-x_min x_text left -x_min right)  xx_text left   1  x_text left xx_text right   fracx_max -xx_max -x_text right -frac12 pi sin left(2 pi fracx_max -xx_max -x_text right right)  xx_text right endcases","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"This is implemented in the following function:","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"FFTLog._c_window","category":"page"},{"location":"components/#FFTLog._c_window","page":"Components","title":"FFTLog._c_window","text":"_c_window(N::AbstractArray, NCut::Int)\n\nThis function returns the smoothing window function as defined in Eq. (C.1) of McEwen et al. (2016).\n\n\n\n\n\n","category":"function"},{"location":"components/","page":"Components","title":"Components","text":"x = Array(0:1000)\ny = FFTLog._c_window(x, 250)\nplot(x,y, ylabel = L\"W(x)\", xlabel = L\"x\")","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"This filtering is found to reduce the ringing effects.","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"After the decomposition, the integral reads","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"F(y)=frac1N sum_m=-N  2^N  2 c_m x_0^nu int_0^infty fracd xx left(fracx_qx_0right)^nu+i eta_m j_ell(x y)","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"Each term is now analytically solvable. ","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"beginaligned\nF(y) =frac1N y^nu sum_m=-N  2^N  2 c_m x_0^-i eta_m y^-i eta_m int_0^infty fracd xx x^nu+i eta_m j_ell(x) = \n=fracsqrtpi4 N y^nu sum_m=-N  2^N  2 c_m x_0^-i eta_m y^-i eta_m g_ellleft(nu+i eta_mright)\nendaligned","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"where g(z) is given by","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"g_ell(z)=2^z fracGammaleft(fracell+z2right)Gammaleft(frac3+ell-z2right) quad-ellRe(z)2","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"giving the range of bias index -ellnu2.","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"Finally, assuming that y is logarithmically sampled with the same linear spacing Delta ln y = Delta x in ln y, we can write the last summation as","category":"page"},{"location":"components/","page":"Components","title":"Components","text":"Fleft(y_pright)=fracsqrtpi4 y_p^nu operatornameIFFTleftc_m^*left(x_0 y_0right)^i eta_m g_ellleft(nu-i eta_mright)right","category":"page"},{"location":"#FFTLog.jl","page":"Home","title":"FFTLog.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FFTLog is a Julia package which performs integrals involving Bessel functions, such as","category":"page"},{"location":"","page":"Home","title":"Home","text":"F(y)=int_0^infty fracd xx f(x) j_ell(x y) ","category":"page"},{"location":"","page":"Home","title":"Home","text":"or Hankel transform","category":"page"},{"location":"","page":"Home","title":"Home","text":"F(y)=int_0^infty mathrmdx f(x) J_n(x y) ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Actually it can evaluate:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Integrals with a Bessel function\nHankel transforms\nIntegrals with a Bessel derivative","category":"page"},{"location":"","page":"Home","title":"Home","text":"In our roadmap, we aim to include also:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Evaluation of C_ℓ's beyond Limber approximation (WIP)\nAutomatic Differentiation\nGPU compatibility\nIntegrals with multiple Bessel functions","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Marco Bonici, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano\nCarmelita Carbone, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the remainder of the documentation, we show how to use FFTLog.jl in details. We also give a quick introduction to the theory behind the FFTLog algorithm.","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please make sure to update tests as appropriate.","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FFTLog.jl is licensed under the MIT \"Expat\" license; see LICENSE for the full license text.","category":"page"}]
}
