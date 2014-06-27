####Operation of analysis:
The run-times discussed here are based on using 256 pixels. If you are using more, you will have to do more on aci by submitting jobs.

Run the 3dft code or the 3dft part of the inverse code to generate the *entire* FT.

You can probably do this on the head node of Odie, hopefully it will only take 5 minutes or less.

Run: ```python lines_to_3d_wave.py ft.gfx ft.txt```
and move the ft.txt file to windows and open in Igor.

Create a new image of 'ft' and add a slider.

Find a spot you want to select (the code will auto-select the mirroring spot)
and record its ```xmin, xman, ymin, ymax, zmin, and zmax``` coordinates (in integers between 0 and npix)
from Igor. Input these numbers into the inverse.f90 code, but make sure to add 1 to each of them because Igor starts at 0 and Fortran starts at 1.

Now run the inverse code (on ACI by submitting a job) while making sure that the filtering is being used, as well as the hanning window.

You should really check that the hanning window looks good before actually running the IFT part,
so run on the Odie head node with the filtered spot and the hanning window and look at the new filtered ft.txt file in Igor.

Mess with the hanning window parameters (```kspotextra = ```... should be all) until you are satisfied.

Now run the full IFT on ACI.

When finished, you may have to copy over the ft.gfx and mgrid.gfx to Odie (probably using rsync because WinSCP is so damn slow!)
(but hopefully you can do the following lines on aci-service-2 since I changed them to not use numpy)
and run: ```python lines_to_3d_wave.py ft.gfx ft.txt && python lines_to_3d_wave.py mgrid.gfx mgrid.txt```

Also run stdev.f90 to run the stdev on mgrid.gfx, but you may have to change 'radius' based on the width of the fringes in Igor.

  'radius' is an integer number representing *half* the fringe width (while AND black).

When stdev is finished (it is pretty short if the radius is around 8 so you can do it on the head node)

run: ```python lines_to_3d_wave.py stdev.gfx stdev.txt``` and copy this to Windows and into Igor. Take a look!

Finally, run:

```python ~/model_analysis/scripts/ift_atom_selection.py Zr50Cu35Al15_t3_final.xyz stdev.gfx```

to generate a 'temp.xyz' file that contains only the atoms you want. Feel free to change the selection criteria at the end of the python file.

The 'temp.xyz' file is usable, so do what you want with it. I do:

```python ~/model_analysis/scripts/ift_cluster_analysis.py Zr50Cu35Al15_t3_final.xyz 3.5 temp.xyz```

where the first model file is the same as the one in the previous command (and what you are using in inverse3dft.f90).

This generates a model called subvpcolored.cif which can be opened with vesta and is a vpcolored model.

That is it! Change as you wish.

Line by line, if you have all the spots indexed that you want to run:

vim inverse3dft.f90 # Change k-spot

make ift && qsub slurm.sh

python lines_to_3d_wave.py mgrid.gfx mgrid.txt # When the ift has finished

python lines_to_3d_wave.py ft.gfx ft.txt # Now copy the files to windows

./stdev

python lines_to_3d_wave.py stdev.gfx stdev.txt # Copy this file too

python ~/model_analysis/scripts/ift_atom_selection.py Zr50Cu35Al15_t3_final.xyz stdev.gfx
