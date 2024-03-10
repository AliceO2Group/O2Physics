#!/bin/bash
script="./run.sh"

$script 211 0.5 pion-50-"$1"
$script 211 0.95 pion-95-"$1"

$script 321 0.5 kaon-50-"$1"
$script 321 0.8 kaon-80-"$1"
$script 321 0.95 kaon-95-"$1"

$script 2212 0.5 proton-50-"$1"
$script 2212 0.8 proton-80-"$1"
$script 2212 0.95 proton-95-"$1"

