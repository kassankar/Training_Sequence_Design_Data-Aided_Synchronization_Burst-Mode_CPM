
# Training Sequence Design for Data-Aided Synchronization of Burst-Mode CPM

A Matlab code for the implementation of the cramer rao bound calculation by __Hosseini & Perrins__  for different types of CPM signals (GMSK-RECT-RC....). The details of the implementation can be obtained from this two papers: _"Training Sequence Design for Data-Aided Synchronization of Burst-Mode CPM"_ and _"The Cramer-Rao Bound for Training Sequence Design for Burst-Mode CPM"_


# How to cite this Work ?
If you use this package for your own work, please consider citing it with [this piece of BibTeX](DATAAIDED_SSB.bib):

```bibtex
@INPROCEEDINGS{9232344,
  author={Kassan, Karim and Farès, Haïfa and Glattli, D. Christian and Louët, Yves},
  booktitle={2020 XXXIIIrd General Assembly and Scientific Symposium of the International Union of Radio Science}, 
  title={Preamble Design for Data-Aided Synchronization of Single Side Band Continuous Phase Modulation}, 
  year={2020},
  volume={},
  number={},
  pages={1-4},
  doi={10.23919/URSIGASS49373.2020.9232344}
  }
```
# How to run the code ?
1. Make sure that you have a compatible version of matlab (this code was tested on matlab 2018b)
2. Download (clone) the files from the repository.
3. Open the file called _CPM_CRB.m_
4. Select the section called __Pulse shape & Variable ini__
5. Select the type of pulses by changing the variable `pulse` number
	* `1` is for Lorentzian
	* `2` is for GMSK
	* `3` is for Raised Cosine
	* `4` is for Rectangular
6. Change the pulse length by changing the variable `L`
7. Select the sampling frequency by changing the variable `os` ( default is `2`)
8. Select M-ary by changing the `M` ( e.g M=2 for Binary)
9. Select the modulation index by changing the `h`.
10. Select the data added preamble length `L0`.
12. Select the section called __Modulated data__
    * select the required SNR values from `Snr`.
## License
© 2020-2021 Karim Kassan

