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
it can evaluate:

- Integrals with a Bessel function
- Hankel transforms
- Integrals with a Bessel derivative

In our roadmap, we aim to include also:
- Evaluation of ``C_{ℓ}``'s beyond Limber approximation (WIP)
- Automatic Differentiation
- GPU compatibility
- Integrals with multiple Bessel functions




### Authors

- Marco Bonici, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano
- Carmelita Carbone, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano


## Usage

In the remainder of the documentation, we show how to use FFTLog.jl in details. We also give
a quick introduction to the theory behind the FFTLog algorithm.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### License

FFTLog.jl is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/marcobonici/FFTLog.jl/blob/main/LICENSE) for
the full license text.
