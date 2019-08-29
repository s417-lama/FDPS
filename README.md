# Task Parallel Implementation for FDPS

Only the calculation of short-range interactions is compatible with task parallelism.

The SPH dambreaking example is in sample/c++/dambreaking.
`Makefile` depends on my environment, so please replace the path of massivethreads.

You can switch the method to parallelize by modifying `Makefile`.

To enable OpenMP, comment out
```
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
```

To enable MassiveThreads, comment out
```
CFLAGS += -DPARTICLE_SIMULATOR_TASK_PARALLEL -I$(MYTH_PATH)/include
LIBS = -lmyth -L $(MYTH_PATH)/lib -Wl,-R$(MYTH_PATH)/lib
```

Run
```
make
./sph.out
```

To configure the settings, modify samples/c++/dambreaking/kernel.h.

* `DAM_2D`: 2d simulation for dambreaking
* `LARGE`: large dataset ONLY for 2d
* `REUSE`: whether to reuse the neighbor list across several iterations
* `OUTPUT`: output the positions of particles per `OUTPUT_INTERVAL`
* `DOUBLE`: double precision or single precision

It uses tpswitch.

# FDPS

FDPS is a general-purpose, high-performance library for particle simulations.

The current version is 5.0b. The previous versions are [here](https://github.com/FDPS/FDPS/releases).

We maintain this from subversion-over-github interface.

If you have some questions, please do not hesitate to contact us. Our
e-mail address is fdps-support-@-mail.jmlab.jp (please replace -@- by @).

Tutorial of FDPS is here
[doc/doc_tutorial_cpp_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_cpp_en.pdf?raw=true)
, Specification is here
[doc/doc_specs_cpp_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_cpp_en.pdf?raw=true).
Tutorial and specification of Fortran/C interface to FDPS are
[doc/doc_tutorial_ftn_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_ftn_en.pdf?raw=true)
and
[doc/doc_specs_ftn_en.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_ftn_en.pdf?raw=true),
respectively.
You can also find a two-page handout here
[doc/doc_SC15_handout.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_SC15_handout.pdf?raw=true).


FDPSのチュートリアルは
[doc/doc_tutorial_cpp_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_cpp_ja.pdf?raw=true)
、仕様書は
[doc/doc_specs_cpp_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_cpp_ja.pdf?raw=true)
にあります。
また、FDPSのFortran/C言語インターフェースのチュートリアルは
[doc/doc_tutorial_ftn_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_tutorial_ftn_ja.pdf?raw=true)
、Fortran/C言語インターフェースの仕様書は
[doc/doc_specs_ftn_ja.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_specs_ftn_ja.pdf?raw=true)
にあります。
2ページでFDPSがわかるハンドアウトはこちらです
[doc/doc_SC15_handout.pdf](https://github.com/FDPS/FDPS/blob/master/doc/doc_SC15_handout.pdf?raw=true)。

ご質問などお問い合わせはfdps-support-@-mail.jmlab.jpにお願いいたします (-@-を@に置換して下さい)。
