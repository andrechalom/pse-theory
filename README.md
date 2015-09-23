# pse - Parameter space exploration

This repository contains the source files used to generate three documents:

* My Master's dissertation, presented at the [Universidade de São Paulo - USP] (http://www.usp.br) (pse.pdf);
* The [Parameter space exploration of ecological models] (http://arxiv.org/abs/1210.6278) paper on arXiv (arxiv1.pdf);
* The [Uncertainty analysis and composite hypothesis under the likelihood paradigm] (http://arxiv.org/abs/1508.03354) paper on arXiv (arxiv2.pdf).

The full text of the dissertation is presented here as the [final.pdf] (http://github.com/andrechalom/pse-theory/raw/master/final.pdf) file. The arxiv1 paper consists on a translation to English of the third chapter and selected sections from the first chapter present in the dissertation, while the arxiv2 paper corresponds to the second chapter.

## Compiling

To generate the pdf version of these documents, download this source code and run the appropriate command:

```
make pse.pdf
make arxiv1.pdf
make arxiv2.pdf
```

This requires: 
* R version 3.0.0 or above, with packages pse, sensitivity, xtable, msm
* A complete Latex environment 
* GNU Make
* Bash

## Other notes

The software `pse` can be found on [this repository] (http://github.com/andrechalom/pse).

This work was supported by a [CAPES] (http://www.capes.gov.br) scholarship (2012-2014).
