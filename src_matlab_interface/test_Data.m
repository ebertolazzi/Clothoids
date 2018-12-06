% Auxiliary values
LU_automatic = 0;
LU_amodio    = 1;
LU_block     = 2;
LU_arceco    = 3;
LU_superLU   = 4;

CS_use_LU           = 0;
CS_use_LUPQ         = 1;
CS_use_QR           = 2;
CS_use_SVD          = 3;
CS_use_LSS          = 4;
CS_use_LSY          = 5;
CS_use_MINIMIZATION = 6;

U_QUADRATIC       = 0;
U_QUADRATIC2      = 1;
U_CUBIC           = 2;
U_QUARTIC         = 3;
U_PARABOLA        = 4;
U_LOGARITHMIC     = 5;
U_COS_LOGARITHMIC = 6;
U_TAN2            = 7;
U_HYPERBOLIC      = 8;
U_BIPOWER         = 9;

PENALTY_REGULAR   = 0;
PENALTY_SMOOTH    = 1;
PENALTY_PIECEWISE = 2;

BARRIER_LOG       = 0;
BARRIER_LOG_EXP   = 1;
BARRIER_LOG0      = 2;

PENALTY2D_RHOMB      = 0;
PENALTY2D_ELLIPSE    = 1;
PENALTY2D_FOURCLOVER = 2;

GREATER_THAN = 0;
LESS_THAN    = 1;
INTERVAL     = 2;
POWER        = 3;
BIPOWER      = 4;

SET  = true;
FREE = false;

data.LU_method             = LU_automatic;
data.ControlSolutionMethod = CS_use_LU;

data.Doctor = false;
data.RedirectStreamToString = false;
data.InfoLevel = 4;
data.OutputSplines = {'s'};

data.Solver.max_iter       = 300;
data.Solver.tolerance      = 1e-09;
data.Solver.min_step       = 0.001;
data.Solver.reduce_factor  = 0.5;
data.Solver.augment_factor = 1.5;
data.Solver.few_iteration  = 8;

data.initial_x = SET;
data.initial_v = SET;
data.final_v   = SET;

data.Controls.FControl.type      = 'U_COS_LOGARITHMIC';
data.Controls.FControl.epsilon   = 0.001;
data.Controls.FControl.tolerance = 0.001;

data.Mesh.s0                 = 0;
data.Mesh.segments{1}.n      = int32(100);
data.Mesh.segments{1}.length = 1;

print_recursive(data);
data
test_gc(data)
