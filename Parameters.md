The input parameters and command to reproduce the dataset of Rc.XDoc dual binding behavior using `MonteCarlo_DBM` are listed as follow,

- Constant pulling speed simulations
```
kT = 4.14 pN*nm
k  = 91   pN/nm
spacer length = 174 nm
```

```
run_simulation(p_w=0.2,num=1000,time=200,step=0.01,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=1.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=40,step=0.001,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=10.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=4.,step=0.0008,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=100.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=2.,step=0.0004,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=400.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=1.,step=0.0002,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=1600.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=0.500,step=0.0001,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=6400.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=0.05,step=0.00001,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=100000.,k=91.,savecurves=False)
run_simulation(p_w=0.2,num=1000,time=0.005,step=0.000001,noiseval=0.,spacer=linker,cpl=xdoc,ldr=load,fingerprint=xMod,forceramp=False,speed=1000000.,k=91.,savecurves=False)
```


- Force clamp
```
run_force_clamp(p_w=0.2, num_sim=1000, cpl=xdoc,fingerprint=xMod, num_force=50, force_from=100, force_interval=10)
```

Information about the input parameters could be found in `##	Simulation command` and `##	Input Parameters` in the py script 
[MonteCarlo_DBM.py](./MonteCarlo_DBM.py)
