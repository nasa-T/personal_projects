"use strict"

var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
var trail = document.getElementById("trail");
var ctx2 = trail.getContext("2d");
const PI = Math.PI;
const dt = 0.01;
const g = 981;
const TIME_SLOWINATOR = 1;
const SHAKE_SIZE = 100;
