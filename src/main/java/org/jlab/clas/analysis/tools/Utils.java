/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis.tools;

import org.jlab.clas.analysis.Particle;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.reactions.DecayKinematics;
import org.jlab.geom.prim.Point3D;

/**
 *
 * @author veronique
 */
public class Utils {
    private double _thetah; //helicity polar angle
    private double _phih; //helicity azimuthal angle
    private double[] _moment; //legendre polynomials moments
    /**
     * 
     * @param parent
     * @param part
     * @return particle is CM frame of parent
     */
    public Particle getCMParticle(Particle parent, Particle part) {
        Particle cmpart = null;
        return cmpart;
    }
    /**
     * 
     * @param event
     * @param part
     * @return particle in event CM frame
     */
    public Particle getCMParticle(Event event, Particle part) {
        Particle cmpart = null;
        return cmpart;
    }

    /**
     * @return the _thetah
     */
    public double getThetah() {
        return _thetah;
    }

    /**
     * @param _thetah the _thetah to set
     */
    public void setThetah(double _thetah) {
        this._thetah = _thetah;
    }

    /**
     * @return the _phih
     */
    public double getPhih() {
        return _phih;
    }

    /**
     * @param _phih the _phih to set
     */
    public void setPhih(double _phih) {
        this._phih = _phih;
    }

    /**
     * @return the _moments
     */
    public double[] getMoments() {
        return _moment;
    }

    /**
     * sets _moment the _moments 
     */
    public void setMoments() {
        int n = 5;
        this._moment = new double[n+1];
        double cq = Math.cos(this._thetah);
        //Bonnet’s recursion formula (n+1)P_{n+1}(x)=(2n+1)xP_{n}(x)-nP_{n-1}(x)
        this._moment[0] = 1;
        this._moment[1] = cq;
        for(int i = 2; i<n; i++) {
            this._moment[i+1] = (1.0/(i+1))*((2*i+1)*cq*this._moment[i]-i*this._moment[i-1]);
        }
    }
    
    public Point3D boost2CMFrm(org.jlab.clas.physics.Particle P, org.jlab.clas.physics.Particle B, LorentzVector Bb) {
        Vector3 boost = Bb.boostVector();
        Vector3 boostm = new Vector3(-boost.x(), -boost.y(), -boost.z());

        LorentzVector  BP1 = new LorentzVector();
        BP1.copy(new LorentzVector(P.px(), P.py(), P.pz(), P.e())); 
        BP1.boost(boostm);
      
        Vector3 partBoosted1 = new Vector3();
        partBoosted1.copy(BP1.vect());
        Vector3 partInFrame1 = DecayKinematics.vectorToFrame(B.vector().vect(), partBoosted1);
        
        return new Point3D(partInFrame1.x(), partInFrame1.y(), partInFrame1.z());
    }
    
    public void boost2CMFrm(org.jlab.clas.physics.Particle P, double beamE, double M) {
        
        org.jlab.clas.physics.Particle B= new org.jlab.clas.physics.Particle(11, 0,0,beamE+M, 0,0,0);
        LorentzVector  Bb = new LorentzVector();
        Bb.copy(B.vector());
        Bb.setPxPyPzE(0.,0., beamE, beamE+M);
        this.boost2CMFrm(P,B, Bb);
        
    }
}
