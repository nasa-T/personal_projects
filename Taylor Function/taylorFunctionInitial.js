var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
var canvas2 = document.getElementById("canvas2");
var ctx2 = canvas2.getContext("2d");
//Increasing DEGREE decreases accuracy the error compunds
const PI = Math.PI;
const DOMAIN_MIN = -4 * PI;
const DOMAIN_MAX = 4 * PI;
const N_PLOTS = 200; //Total among all graphs
const DOMAIN_LENGTH = DOMAIN_MAX - DOMAIN_MIN;
const RADIUS = 0.5;
const SCALE = canvas.width / (DOMAIN_LENGTH);
const dt = 0.01;
const DEGREE_CHANGE = 2;
const MAX_DEGREE = 20; //Don't go past 20.
const SHIFT = 0; //Where to model Taylor Series around

var DEGREE = 0; //How many terms in Taylor Series? Starts code DEGREE_CHANGE after DEGREE
var x = DOMAIN_MIN;
