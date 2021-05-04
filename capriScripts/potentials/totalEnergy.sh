#!/bin/bash

model=$1 # model PDB
rchain=$2 # receptor chain
lchain=$3 # ligand chain

scriptsDir=$HOME/scripts/potentials

pisaEnergy=`$scriptsDir'/pisa/pisaEnergy' $model $rchain $lchain $scriptsDir'/pisa/pisa.params'`
pieEnergy=`$scriptsDir'/pie/bin/pie_score' $model $rchain $lchain $scriptsDir'/pie/bin/pie.params'`
hbEnergy=`$scriptsDir'/hbnn/hbnn' $model $rchain $lchain $scriptsDir'/hbnn/hbnet.params'`

c3Energy=`echo $pieEnergy' '$pisaEnergy | awk '{print -0.8*$pieEnergy + 0.1*$pisaEnergy +$pieEnergy*$pisaEnergy}'`
totalEnergy=`echo $c3Energy' '$hbEnergy | awk '{print $c3Energy + 3.3*$hbEnergy}'`

echo $totalEnergy
