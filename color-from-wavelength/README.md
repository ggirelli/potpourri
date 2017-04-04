color-from-wavelength
===

I needed a function to convert a wavelength into a color in hex format. I started googling and stumbled upon [this page](http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm). There, they present a nice function in Borland Pascal, which in turn is based on the Fortran function that can be found [here](http://www.midnightkite.com/color.html).

Since I needed it for both back-end and front-end of a Python Bottle-based webserver with an Angular2 framework, I re-implemented them in TypeScript and in Python.

Both functions are, as those they are based on, quite straightforward and can be easily customized to get the output in RGB format.

Cheers!

└[∵┌]└[ ∵ ]┘[┐∵]┘