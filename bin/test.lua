-- Dati test N1.

DATA = {

  mesh_files = { "granny_toro.1.neutral.bz2",
                 "granny_toro.2.neutral.bz2",
                 "granny_toro.3.neutral.bz2",
                 "granny_toro.4.neutral.bz2" },

  do_save  = { true, true, true, false },
  integers = { 1, 2, 3, 4, 3, 2, 1 },
  reals    = { 1, 2.1, 3, 4, 3, 2, 1 },

  out_dir = "granny_toro/",

  regions = {
    { "toro", 1 },
    { "granny", 2 },
    { "insulant", 3 }
  },

  mu    = 1.25e-6,
  omega = 50*2*3.14,
  sigma = 4.0e7,

  a  = 3.88,
  b  = 2.57,
  c  = 2.09,

}
