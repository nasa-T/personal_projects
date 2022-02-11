var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
var canvas2 = document.getElementById("canvas2");
var ctx2 = canvas2.getContext("2d");
//Increasing DEGREE decreases accuracy the error compunds
const PI = Math.PI;
const DOMAIN_MIN = -4 * PI;
const DOMAIN_MAX = 4 * PI;
const N_PLOTS = 300; //Total among all graphs
const DOMAIN_LENGTH = DOMAIN_MAX - DOMAIN_MIN;
const RADIUS = 0.5;
const SCALE = canvas.width / (DOMAIN_LENGTH);
const DEGREE_CHANGE = 1;
const MAX_DEGREE = 12; //Don't go past 20.
const RANDOM_NUM_RANGE = 10;

var SHIFT = DOMAIN_MIN; //Where to model Taylor Series around

var DEGREE = 1; //How many terms in Taylor Series? Starts code DEGREE_CHANGE after DEGREE
var dt = 0.1;
