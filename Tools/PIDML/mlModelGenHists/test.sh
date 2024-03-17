#!/bin/bash
script="./run.sh"

$script 211 0.5 pion-plus-50-"$1"
$script 211 0.95 pion-plus-95-"$1"

$script -211 0.5 pion-minus-50-"$1"
$script -211 0.95 pion-minus-95-"$1"

$script 321 0.5 kaon-plus-50-"$1"
$script 321 0.8 kaon-plus-80-"$1"
$script 321 0.95 kaon-plus-95-"$1"

$script -321 0.5 kaon-minus-50-"$1"
$script -321 0.8 kaon-minus-80-"$1"
$script -321 0.95 kaon-minus-95-"$1"

$script 2212 0.5 proton-plus-50-"$1"
$script 2212 0.8 proton-plus-80-"$1"
$script 2212 0.95 proton-plus-95-"$1"

$script -2212 0.5 proton-minus-50-"$1"
$script -2212 0.8 proton-minus-80-"$1"
$script -2212 0.95 proton-minus-95-"$1"

