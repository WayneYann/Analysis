from ase import Atoms
from ase.io.trajectory import Trajectory as traj
from ase.lattice import surface
from ase.constraints import FixAtoms
from ase.optimize import LBFGS
from ase.calculators.vasp import Vasp
from ase.io.vasp import write_vasp
from ase.visualize import view
from ase.io import write

calc = Vasp(xc='PW91', kpts=(1,1,1), nwrite=1, lwave=False, lcharg=False,lvtot=False
            , encut=400, algo='Fast', ismear=1, sigma=0.01, voskown=1, istart=0, nelm=400, nelmdl=-10, ediff=1e-6, ispin=2
            ,nsw=1000, isif=2, ibrion=1, nfree=2, potim=0.2, ediffg=-0.04
            ,npar=1
            ,lreal='Auto')

# Create a c(2x2) surface with 4 layers and 14 Angstrom of vacuum
d=1.5
molecule = Atoms('O2', positions=[(0,0,0), (0,0,d)])

molecule.center(vacuum=14.0)
view(molecule) # View the slab, frozen atoms will be marked with a "X"


#slab.set_calculator(calc)

calc.initialize(molecule)
calc.write_incar(molecule)
calc.write_potcar()
calc.write_kpoints()
write('POSCAR', calc.atoms_sorted)  # this will write a "sorted" POSCAR

