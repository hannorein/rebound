typedef struct {
    double ID;
    double x, y, z;
    double vx, vy, vz;
    double mass, density;
} ReadParticle;

// Function declarations
void transform(ReadParticle *p, double x_translation);
