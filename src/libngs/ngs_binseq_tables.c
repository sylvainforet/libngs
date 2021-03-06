/* Copyright (C) 2010  Sylvain FORET
 *
 * This file is part of libngs.
 *
 * libngs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *
 */

const unsigned char revcomp_table[256] =
{
  [  0] = 255,
  [  1] = 191,
  [  2] = 127,
  [  3] =  63,
  [  4] = 239,
  [  5] = 175,
  [  6] = 111,
  [  7] =  47,
  [  8] = 223,
  [  9] = 159,
  [ 10] =  95,
  [ 11] =  31,
  [ 12] = 207,
  [ 13] = 143,
  [ 14] =  79,
  [ 15] =  15,
  [ 16] = 251,
  [ 17] = 187,
  [ 18] = 123,
  [ 19] =  59,
  [ 20] = 235,
  [ 21] = 171,
  [ 22] = 107,
  [ 23] =  43,
  [ 24] = 219,
  [ 25] = 155,
  [ 26] =  91,
  [ 27] =  27,
  [ 28] = 203,
  [ 29] = 139,
  [ 30] =  75,
  [ 31] =  11,
  [ 32] = 247,
  [ 33] = 183,
  [ 34] = 119,
  [ 35] =  55,
  [ 36] = 231,
  [ 37] = 167,
  [ 38] = 103,
  [ 39] =  39,
  [ 40] = 215,
  [ 41] = 151,
  [ 42] =  87,
  [ 43] =  23,
  [ 44] = 199,
  [ 45] = 135,
  [ 46] =  71,
  [ 47] =   7,
  [ 48] = 243,
  [ 49] = 179,
  [ 50] = 115,
  [ 51] =  51,
  [ 52] = 227,
  [ 53] = 163,
  [ 54] =  99,
  [ 55] =  35,
  [ 56] = 211,
  [ 57] = 147,
  [ 58] =  83,
  [ 59] =  19,
  [ 60] = 195,
  [ 61] = 131,
  [ 62] =  67,
  [ 63] =   3,
  [ 64] = 254,
  [ 65] = 190,
  [ 66] = 126,
  [ 67] =  62,
  [ 68] = 238,
  [ 69] = 174,
  [ 70] = 110,
  [ 71] =  46,
  [ 72] = 222,
  [ 73] = 158,
  [ 74] =  94,
  [ 75] =  30,
  [ 76] = 206,
  [ 77] = 142,
  [ 78] =  78,
  [ 79] =  14,
  [ 80] = 250,
  [ 81] = 186,
  [ 82] = 122,
  [ 83] =  58,
  [ 84] = 234,
  [ 85] = 170,
  [ 86] = 106,
  [ 87] =  42,
  [ 88] = 218,
  [ 89] = 154,
  [ 90] =  90,
  [ 91] =  26,
  [ 92] = 202,
  [ 93] = 138,
  [ 94] =  74,
  [ 95] =  10,
  [ 96] = 246,
  [ 97] = 182,
  [ 98] = 118,
  [ 99] =  54,
  [100] = 230,
  [101] = 166,
  [102] = 102,
  [103] =  38,
  [104] = 214,
  [105] = 150,
  [106] =  86,
  [107] =  22,
  [108] = 198,
  [109] = 134,
  [110] =  70,
  [111] =   6,
  [112] = 242,
  [113] = 178,
  [114] = 114,
  [115] =  50,
  [116] = 226,
  [117] = 162,
  [118] =  98,
  [119] =  34,
  [120] = 210,
  [121] = 146,
  [122] =  82,
  [123] =  18,
  [124] = 194,
  [125] = 130,
  [126] =  66,
  [127] =   2,
  [128] = 253,
  [129] = 189,
  [130] = 125,
  [131] =  61,
  [132] = 237,
  [133] = 173,
  [134] = 109,
  [135] =  45,
  [136] = 221,
  [137] = 157,
  [138] =  93,
  [139] =  29,
  [140] = 205,
  [141] = 141,
  [142] =  77,
  [143] =  13,
  [144] = 249,
  [145] = 185,
  [146] = 121,
  [147] =  57,
  [148] = 233,
  [149] = 169,
  [150] = 105,
  [151] =  41,
  [152] = 217,
  [153] = 153,
  [154] =  89,
  [155] =  25,
  [156] = 201,
  [157] = 137,
  [158] =  73,
  [159] =   9,
  [160] = 245,
  [161] = 181,
  [162] = 117,
  [163] =  53,
  [164] = 229,
  [165] = 165,
  [166] = 101,
  [167] =  37,
  [168] = 213,
  [169] = 149,
  [170] =  85,
  [171] =  21,
  [172] = 197,
  [173] = 133,
  [174] =  69,
  [175] =   5,
  [176] = 241,
  [177] = 177,
  [178] = 113,
  [179] =  49,
  [180] = 225,
  [181] = 161,
  [182] =  97,
  [183] =  33,
  [184] = 209,
  [185] = 145,
  [186] =  81,
  [187] =  17,
  [188] = 193,
  [189] = 129,
  [190] =  65,
  [191] =   1,
  [192] = 252,
  [193] = 188,
  [194] = 124,
  [195] =  60,
  [196] = 236,
  [197] = 172,
  [198] = 108,
  [199] =  44,
  [200] = 220,
  [201] = 156,
  [202] =  92,
  [203] =  28,
  [204] = 204,
  [205] = 140,
  [206] =  76,
  [207] =  12,
  [208] = 248,
  [209] = 184,
  [210] = 120,
  [211] =  56,
  [212] = 232,
  [213] = 168,
  [214] = 104,
  [215] =  40,
  [216] = 216,
  [217] = 152,
  [218] =  88,
  [219] =  24,
  [220] = 200,
  [221] = 136,
  [222] =  72,
  [223] =   8,
  [224] = 244,
  [225] = 180,
  [226] = 116,
  [227] =  52,
  [228] = 228,
  [229] = 164,
  [230] = 100,
  [231] =  36,
  [232] = 212,
  [233] = 148,
  [234] =  84,
  [235] =  20,
  [236] = 196,
  [237] = 132,
  [238] =  68,
  [239] =   4,
  [240] = 240,
  [241] = 176,
  [242] = 112,
  [243] =  48,
  [244] = 224,
  [245] = 160,
  [246] =  96,
  [247] =  32,
  [248] = 208,
  [249] = 144,
  [250] =  80,
  [251] =  16,
  [252] = 192,
  [253] = 128,
  [254] =  64,
  [255] =   0
};

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
