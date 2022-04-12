# FFTLog.jl

FFTLog is a Julia package which performs integrals involving Bessel functions, such as

```math
F(y)=\int_{0}^{\infty} \frac{d x}{x} f(x) j_{\ell}(x y)\, ,
```

or Hankel transform

```math
F(y)=\int_{0}^{\infty} \mathrm{d}x f(x) J_{n}(x y)\, .
```

Actually
it can handle:

- Integrals with a Bessel function
- Integrals with a Bessel derivative
- Hankel transforms

In our roadmap, we aim to include also:
- Automatic Differentiation
- GPU compatibility
- Integrals with multiple Bessel functions




### Authors

- Marco Bonici, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano


## Usage

In the remainder of the documentation, we explain the theory behind the FFTLog algorithm.
However, using FFTLog.jl is quite easy.

1. Instantiate an object `HankelPlan`
```julia
HankelTest = FFTLog.HankelPlan(x = k)
```
2. Perform some precomputations
```julia
prepare_Hankel!(HankelTest, Ell)
```
3. Compute the Hankel transform
```julia
Fy = evaluate_Hankel(HankelTest, Pk)
```
4. If needed, the array `y` (the counterpart of the array `r`) can be obtained with
```julia
y = get_y(HankelTest)
```




## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### License

FFTLog.jl is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/marcobonici/FFTLog.jl/blob/main/LICENSE) for
the full license text.
