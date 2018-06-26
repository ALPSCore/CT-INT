# List of tutorials
- [Tutorial1](#tutorial1)
- [Tutorial2](#tutorial2)

## <a ref="tutorial1">Tutorial1</a>
In this tutorial, we solve the three-orbital impurity model described in our preprint.

You will see how to run the solver and process the output data in job.sh.

```
```

At the end of the simulation, the results are written into the HDF5 file "input.out.h5".
A HDF5 file can be read easily e.g. by using the h5py library in Python.
By running "read_G.py", you can the QMC data into text files.
Or, you can plot the self-energy by running "plot.py".
![](tutorial1/Sigma-Re.png)
![](tutorial1/Sigma-Im.png)


## <a ref="tutorial2">Tutorial2</a>
In tutorial 2, we solve exactly the same model as in tutorial 1 by combining spin-up and spin-down components into a single flavor.
This sample will be usuful when you want to solve a model with spin-orbit interaction.
