class plasMDSolver : public Molecule {
public:
    plasMDSolver();
    void set_beam_profile();
    void Solve();
private:
    void timestep(double dt);
};
