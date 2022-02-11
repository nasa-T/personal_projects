var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");

const EARTH_ORBIT = 1.5e11;
const EARTH_MASS = 6e24;
const SUN_MASS = 2e30;
const G = 6.67 * Math.pow(10, -11);
const PI = Math.PI;
const N_G = 40;
const N_B = 1;
const GRID_WIDTH = 5e11;
const GRID_HEIGHT = 5e11;

const N_TAIL = 100;
const dt = 24 * 60 * 60;
const YEAR = 365.24 * 24 * 60 * 60;
const SCALE_TIME = 20 / YEAR;
var SCALE = canvas.width / (GRID_WIDTH);
var PAN_X = 0;
var PAN_Y = 0;
