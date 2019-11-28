#! /bin/bash
########################################################
GMXDIR=/pubhome/yzhou/Software/silcsbio/gromacs20183/bin
numsys=2
frag_density=0.0006
########################################################

cp data/*.dat ./
./data/make_pdb_charmm_compatible charmm36.ff rec.pdb rec_4pdb2gmx.pdb
grep -vE 'TIP|HOH' rec_4pdb2gmx.pdb > tmp.pdb
grep -E 'TIP|HOH'  rec_4pdb2gmx.pdb >> tmp.pdb
mv tmp.pdb rec_4pdb2gmx.pdb

string=`./data/generate_termini_caps_strings rec_4pdb2gmx.pdb`
echo $string | $GMXDIR/gmx pdb2gmx -f rec_4pdb2gmx.pdb -o rec_gmx.pdb -p rec_gmx.top -ff charmm36 -water tip3p -ter -merge all &> gmx.log
if [[ ! -f "rec_gmx.pdb" ]];
  then
    echo "Gromacs utility PDB2GMX did not generate the pdb/top !"
    exit
fi

echo "q" | $GMXDIR/gmx make_ndx -f rec_gmx.pdb -o selection.ndx
echo "3" | $GMXDIR/gmx genrestr -f rec_gmx.pdb -fc 50.208 50.208 50.208 -o posre_protein_ca.itp -n selection.ndx

sed '/; Include water topology/i \; Include C-alpha positional restraints file\n#ifdef POSRE_CA\n\#include \"posre_protein_ca.itp\"\n\#endif\n' rec_gmx.top > tmp.top
mv tmp.top rec_gmx.top
sed '/; Include water topology/i \; Include silcs molecules\n\#include \"./charmm36.ff/mol/clbx.itp\"\n' rec_gmx.top > tmp.top
mv tmp.top rec_gmx.top

i=1
  while [[ "$i" -le ${numsys} ]];
  do
     $GMXDIR/gmx editconf -f rec_gmx.pdb -o rec_gmx_box_${i}.pdb -d 1.5 &>> gmx.log

     xtala=`grep -e CRYST1 rec_gmx_box_${i}.pdb | awk '{printf "%1.3f", $2}'`
     xtalb=`grep -e CRYST1 rec_gmx_box_${i}.pdb | awk '{printf "%1.3f", $3}'`
     xtalc=`grep -e CRYST1 rec_gmx_box_${i}.pdb | awk '{printf "%1.3f", $4}'`
     num_per_box=`echo 1 | awk -v a=$xtala -v b=$xtalb -v c=$xtalc -v den=$frag_density '{printf "%d",a*b*c*den}'`
     echo "${num_per_box} molecules per box will be added before deletion"

     cp rec_gmx.top rec_solv.${i}.top
     ${GMXDIR}/gmx insert-molecules -seed ${i}31415 -f rec_gmx_box_${i}.pdb -ci charmm36.ff/mol/clbx.pdb -o temp${i}.pdb -nmol $num_per_box &>> gmx.log
     marker=`grep "^ATOM" charmm36.ff/mol/clbx.pdb | awk '{print $3}' | uniq | head -1`
     nfrag=`grep "$(echo clbx | tr '[:lower:]' '[:upper:]' | cut -c1-3)" temp${i}.pdb | grep ${marker} | wc -l`
     echo "$(echo clbx | tr '[:lower:]' '[:upper:]')     ${nfrag}" >> rec_solv.${i}.top

     # add water
     ${GMXDIR}/gmx solvate -cp temp${i}.pdb -cs -o rec_solv.${i}.pdb -p rec_solv.${i}.top &>> gmx.log

     # add ions
     ${GMXDIR}/gmx grompp -f data/ions.mdp -c rec_solv.${i}.pdb -p rec_solv.${i}.top -o ions.${i}.tpr  &>> gmx.log
     ${GMXDIR}/gmx genion -s ions.${i}.tpr -o rec_solv_ions.${i}.pdb -p rec_solv.${i}.top -pname NA -nname CL -neutral

     if [[ ! -d MD-$i ]]; then mkdir MD-$i; fi
     cp -r rec_solv_ions.$i.pdb rec_solv.$i.top *.itp charmm36.ff MD-$i/
     echo "#\$ -S /bin/bash" > MD-$i/job_MD.cmd
     echo "#\$ -cwd" >> MD-$i/job_MD.cmd
     echo "#\$ -V" >> MD-$i/job_MD.cmd
     echo "#\$ -o MD.out" >> MD-$i/job_MD.cmd
     echo "#\$ -e MD.err" >> MD-$i/job_MD.cmd
     echo "#\$ -j y" >> MD-$i/job_MD.cmd
     echo "#\$ -N MD_${i}" >> MD-$i/job_MD.cmd
     echo "#\$ -pe cuda 1" >> MD-$i/job_MD.cmd
     echo "#\$ -l gpu=1" >> MD-$i/job_MD.cmd
     echo "export OMP_NUM_THREAD=2" >> MD-$i/job_MD.cmd
     echo "source /usr/bin/startcuda.sh" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx grompp -f ../data/emin.mdp -r rec_solv_ions.$i.pdb -c rec_solv_ions.$i.pdb -p rec_solv.$i.top -o em.tpr" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx mdrun -v -deffnm em" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx grompp -f equil.mdp  -c em.gro -r em.gro -p rec_solv.$i.top -o nvt.tpr" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx mdrun -v -deffnm nvt" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx grompp -f prod.mdp -r nvt.gro -c nvt.gro -p rec_solv.$i.top -o nptMD.tpr" >> MD-$i/job_MD.cmd
     echo "$GMXDIR/gmx mdrun -deffnm nptMD" >> MD-$i/job_MD.cmd
     echo "source /usr/bin/end_cuda.sh" >> MD-$i/job_MD.cmd

     i=$((i+1))
  done

rm -f temp*.pdb
