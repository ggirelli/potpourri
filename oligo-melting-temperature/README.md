oligo-melting-temperature
===

The `oligo_tm_calc.py` script, implemented in Python, allows to calculate the melting temperature of a DNA duplex, provided the sequence of one of the two strands.

The calculation is based on the N-N thermodynamic values presented in Allawi&Santalucia[1]. The melting temperature calculation is based on Santalucia, 1998[2]. Sodium and cagnesium concentration correction is based on the work of Owczarzy et al[3,4].

By specifying the `-t` option, it is also possible to use the N-N thermodynamic values presented in Freier et al[5] for RNA.

By using the `-F` option and providing the path to a file with one oligo sequence per line, the melting temperature is automatically calculated for every sequence in the file.

### Help page

```
usage: oligo_tm_calc.py [-h] [-o oligo_conc] [-n na_conc] [-m mg_conc]
                        [-t {allawi,freier}] [-C] [-F] [-v]
                        seq

Calculate melting temeprature of a DNA duplex at provided [oligo], [Na+],
[Mg2+]. Either provide an oligo sequence or a file with one oligo per line
(and use -F option). References: [1] Freier et al, PNAS(83), 1986; [2] Allawi
& Santalucia, Biochemistry(36), 1997; [3] SantaLucia, PNAS(95), 1998; [4]
Owczarzy et al, Biochemistry(43), 2004; [5] Owczarzy et al, Biochemistry(47),
2008.

positional arguments:
  seq                   DNA duplex sequence (one strand only) or path to file
                        containing one sequence per line (use with -F).

optional arguments:
  -h, --help            show this help message and exit
  -o oligo_conc, --oconc oligo_conc
                        Oligonucleotide concentration [M]. Default: 0.25e-6 M
  -n na_conc, --naconc na_conc
                        Na+ concentration [M]. Default: 50e-3 M
  -m mg_conc, --mgconc mg_conc
                        Mg2+ concentration [M]. Default: 0 M
  -t {allawi,freier}, --thermo {allawi,freier}
                        Thermodynamic table to use in the calculations.
                        Possible values: allawi (based on ref.2, default) or
                        freier (based on ref.1).
  -C, --celsius         Output temperature in Celsius degrees. Default: Kelvin
  -F, --usefile         Use when a file path is provided instead of a single
                        sequence.
  -v, --verbose         Verbose output.

```

### References

1. Allawi & Santalucia, Biochemistry(36), 1997
2. SantaLucia, PNAS(95), 1998
3. Owczarzy et al, Biochemistry(43), 2004
4. Owczarzy et al, Biochemistry(47), 2008
5. Freier et al, PNAS(83), 1986