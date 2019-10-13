#ifndef __CHARGWD_PARTICLE_TRACER_
#define __CHARGWD_PARTICLE_TRACER_

#include <vector>
#include <array>
#include <string>


class ChargedParticleTracer {
public:
    void ReserveMemory(std::size_t numberOfElements) {
        charges.reserve(numberOfElements);
        masses.reserve(numberOfElements);
        positions.reserve(numberOfElements);
        velocities.reserve(numberOfElements);
        momentums.reserve(numberOfElements);
        forces.reserve(numberOfElements);
    }

    void AddParticle(const double charge,
                     const double mass,
                     const std::array<double, 3>& position,
                     const std::array<double, 3>& velocity,
                     const std::array<double, 3>& force) {
        double velocitySquared = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];
        double gamma = 1.0/std::sqrt(1 - velocitySquared/c_squared);
        std::array<double, 3> momentum{gamma*mass*velocity[0],
                                       gamma*mass*velocity[1],
                                       gamma*mass*velocity[2]};
        charges.push_back(charge);
        masses.push_back(mass);
        positions.push_back(position);
        velocities.push_back(velocity);
        momentums.push_back(momentum);
        forces.push_back(force);
    }

    void UpdateParticlesMomentumsAndVelocities(const double dt) {
        for(std::size_t i = 0; i < positions.size(); ++i) {
            std::array<double, 3>& p = momentums[i];
            std::array<double, 3>& f = forces[i];
            p[0] += f[0]*dt;
            p[1] += f[1]*dt;
            p[2] += f[2]*dt;

            double mSquared = masses[i]*masses[i];
            double pSquared = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
            std::array<double, 3>& v = velocities[i];
            double mRel = std::sqrt(mSquared + pSquared/c_squared);
            //mRel = masses[i];
            v[0] = p[0]/mRel;
            v[1] = p[1]/mRel;
            v[2] = p[2]/mRel;
        }
    }

    void UpdateParticlesPositions(const double dt) {
        for(std::size_t i = 0; i < positions.size(); ++i) {
            std::array<double, 3>& r = positions[i];
            std::array<double, 3>& v = velocities[i];
            r[0] += v[0]*dt;
            r[1] += v[1]*dt;
            r[2] += v[2]*dt;
        }
    }

    void UpdateParticles(const double dt) {
        UpdateParticlesMomentumsAndVelocities(dt);
        UpdateParticlesPositions(dt);
    }

    void SaveData(std::string folder, int timeInd, double time,
                  std::string filename_masscharge="masscharge",
                  std::string filename_position_velocity_momentum="posvelmom"
                  ) {
        std::string fileSuffix = std::string("_")                                  \
                                 + "_i=" + boost::lexical_cast<std::string>(timeInd)           \
                                 + "_t=" + boost::lexical_cast<std::string>(time);

        std::string fileName = folder + "/" + filename_masscharge + fileSuffix;
        std::ofstream fileOut(fileName.c_str(), std::ios::out | std::ios::binary);
        assert(fileOut.is_open());

        fileOut.write((char*)(masses.data()), masses.size() * sizeof(double));
        fileOut.write((char*)(charges.data()), charges.size() * sizeof(double));
        fileOut.close();

        std::string fileName_p = folder + "/" + filename_position_velocity_momentum + fileSuffix;
        std::ofstream fileOut_p(fileName_p.c_str(), std::ios::out | std::ios::binary);
        assert(fileOut_p.is_open());

        for (int i = 0; i < positions.size(); ++i) {
            fileOut_p.write((char*)(&(positions[i][0])), 3 * sizeof(double));
        }
        for (int i = 0; i < velocities.size(); ++i) {
            fileOut_p.write((char*)(&(velocities[i][0])), 3 * sizeof(double));
        }
        for (int i = 0; i < momentums.size(); ++i) {
            fileOut_p.write((char*)(&(momentums[i][0])), 3 * sizeof(double));
        }
        fileOut_p.close();

    }


    auto& GetMasses() {return masses;}
    auto& GetCharges() {return charges;}
    auto& GetPositions() {return positions;}
    auto& GetVelocities() {return velocities;}
    auto& GetMomentums() {return momentums;}
    auto& GetForces() {return forces;}

private:
    double c_squared = PhysicalConstants_SI::speedOfLight * PhysicalConstants_SI::speedOfLight;
    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<std::array<double, 3>> positions;
    std::vector<std::array<double, 3>> velocities;
    std::vector<std::array<double, 3>> momentums;
    std::vector<std::array<double, 3>> forces;

};

#endif // __CHARGWD_PARTICLE_TRACER_
