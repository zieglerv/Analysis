/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import org.jlab.clas.swimtools.Swim;
import cnuphys.swim.SwimTrajectory;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
/**
 *
 * @author ziegler
 */
public class Particle {

    public boolean isEBParticle;
    private int[] _ndau;

    public Particle(Particle part1) {
        this._E=part1._E;
        this._Idx = part1._Idx;
        this._beta = part1._beta;
        this._calcBeta = part1._calcBeta;
        this._charge=part1._charge;
        this._det = part1._det;
        this._daughters = part1._daughters;
        this._mass = part1._mass;
        this._massConstrE = part1._massConstrE;
        this._mass = part1._mass;
        this._umass = part1._umass;
        this._uncormass = part1._uncormass;
        this._p = part1._p;
        this._pid = part1._pid;
        this._pt = part1._pt;
        this._px = part1._px;
        this._py = part1._py;
        this._pz = part1._pz;
        this._upx = part1._upx;
        this._upy = part1._upy;
        this._upz = part1._upz;
        this._r = part1._r;
        this._recmass = part1._recmass;
        this._status = part1._status;
        this._uE = part1._uE;
        this._massConstrE = part1._massConstrE;
        this._ovx = part1._ovx;
        this._ovy = part1._ovy;
        this._ovz = part1._ovz;
        this._vx = part1._vx;
        this._vy = part1._vy;
        this._vz = part1._vz;
    }
    public double massConstrE() {
        double e = this.getMassConstrE();
        if(Constants.particleMap.containsKey(this.getPid())) {
            double cmass = Constants.particleMap.get(Math.abs(this.getPid())).mass();
            e = Math.sqrt(this.getP()*this.getP()+cmass*cmass);
        }
        return e;
    }
    /**
     * @return the recParticle
     */
    public org.jlab.clas.physics.Particle getRecParticle() {
        return recParticle;
    }

    /**
     * @param recParticle the recParticle to set
     */
    public void setRecParticle(org.jlab.clas.physics.Particle recParticle) {
        this.recParticle = recParticle;
    }
    
    public boolean isUsed = false;
    private double _mass ;
    private double _recmass ;
    private double _E;
    private double _massConstrE;
    private double _px;
    private double _py;
    private double _pz;
    private double _pxcm;
    private double _pycm;
    private double _pzcm;
    private double _vx;
    private double _vy;
    private double _vz;
    private double _r; // distance of closest approach between 2 track helices
    private double _p;
    private double _pt;
    private double _uncormass;
    private double _umass ;
    private double _uE;
    private double _upx;
    private double _upy;
    private double _upz;
    private double _ovx;
    private double _ovy;
    private double _ovz;
    
    private double _beta;
    private double _calcBeta;
    private int _charge;
    private int _status;
    private int _Idx=-1;
    private int _pid;
    public int vIndex=-1; 
    
    private int _det = -1;
    
    private List<Particle> _daughters = new ArrayList<Particle>();
    private org.jlab.clas.physics.Particle recParticle ;
    
    public Particle() {
        _daughters = new ArrayList<Particle>();
    }

    /**
     * @return the _ndau
     */
    public int[] getNdau() {
        return _ndau;
    }

    /**
     * @param _ndau the _ndau to set
     */
    public void setNdau(int[] _ndau) {
        this._ndau = _ndau;
    }
    
    private void setEBParticle(double beta, int pid, double upx, double upy, double upz,
            double vx, double vy, double vz, int charge) {
        _daughters = new ArrayList<Particle>();
        this._beta = beta;
        this._ovx = vx;
        this._ovy = vy;
        this._ovz = vz;
        this._upx = upx;
        this._upy = upy;
        this._upz = upz;
        this._pid = pid;
        this._p = Math.sqrt(upx*upx+upy*upy+upz*upz);
        this._pt = Math.sqrt(upx*upx+upy*upy);
        double cmass = Constants.particleMap.get(Math.abs(pid)).mass();
        double mass = Math.sqrt((this._p/beta)*(this._p/beta)-this._p*this._p);
        double calcBeta = this._p/Math.sqrt(this._p*this._p+cmass*cmass);
        this._charge = charge;
        this._mass = cmass;
        this._umass = mass;
        this._calcBeta = calcBeta;
        double E = Math.sqrt(this._p*this._p+mass*mass);
        double cE = Math.sqrt(this._p*this._p+cmass*cmass);
        this._uE =E;
        this._E = cE;  //always use final states mass-constrained energy
        this._massConstrE = cE;
        this.isEBParticle = true; 
    }
    private void setPhoton(double energy, double beta, double upx, double upy, double upz,
            double vx, double vy, double vz, int charge) {
        _daughters = new ArrayList<Particle>();
        this._beta =beta;
        this._calcBeta =1;
        this._ovx = vx;
        this._ovy = vy;
        this._ovz = vz;
        this._upx = upx;
        this._upy = upy;
        this._upz = upz;
        this._pid = 22;
        this._p = Math.sqrt(upx*upx+upy*upy+upz*upz);
        this._pt = Math.sqrt(upx*upx+upy*upy);
        double cmass = Constants.particleMap.get(Math.abs(22)).mass();
        
        this._charge = charge;
        this._mass = cmass;
        double m2 = energy*energy-this._p*this._p;
        if(m2<0) m2 = 0;
        this._umass = Math.sqrt(m2);
        this._E =energy;
        this._massConstrE = Math.sqrt(cmass*cmass+this._p*this._p);
        this.isEBParticle = true;
    }
    
    
    
    
    public Particle(double beta, int pid, double upx, double upy, double upz,
            double vx, double vy, double vz, int charge) {
        this.setEBParticle(beta, pid, upx, upy, upz, vx, vy, vz, charge);
    }    

    public Particle(double energy, double beta, double upx, double upy, double upz,
            double vx, double vy, double vz, int charge) {
        this.setPhoton(energy, beta, upx, upy, upz, vx, vy, vz, charge);
    }
    /**
     * @return the _umass
     */
    public double getUmass() {
        return _umass;
    }

    /**
     * @param _umass the _umass to set
     */
    public void setUmass(double _umass) {
        this._umass = _umass;
    }

    /**
     * @return the _uE
     */
    public double getuE() {
        return _uE;
    }

    /**
     * @param _uE the _uE to set
     */
    public void setuE(double _uE) {
        this._uE = _uE;
    }

    /**
     * @return the _upx
     */
    public double getUpx() {
        return _upx;
    }

    /**
     * @param _upx the _upx to set
     */
    public void setUpx(double _upx) {
        this._upx = _upx;
    }

    /**
     * @return the _upy
     */
    public double getUpy() {
        return _upy;
    }

    /**
     * @param _upy the _upy to set
     */
    public void setUpy(double _upy) {
        this._upy = _upy;
    }

    /**
     * @return the _upz
     */
    public double getUpz() {
        return _upz;
    }

    /**
     * @param _upz the _upz to set
     */
    public void setUpz(double _upz) {
        this._upz = _upz;
    }

    /**
     * @return the _uvx
     */
    public double getOvx() {
        return _ovx;
    }

    /**
     * @param _uvx the _uvx to set
     */
    public void setOvx(double _uvx) {
        this._ovx = _uvx;
    }

    /**
     * @return the _uvy
     */
    public double getOvy() {
        return _ovy;
    }

    /**
     * @param _uvy the _uvy to set
     */
    public void setOvy(double _uvy) {
        this._ovy = _uvy;
    }

    /**
     * @return the _uvz
     */
    public double getOvz() {
        return _ovz;
    }

    /**
     * @param _uvz the _uvz to set
     */
    public void setOvz(double _uvz) {
        this._ovz = _uvz;
    }

    /**
     * @return the _det
     */
    public int getDet() {
        return _det;
    }

    /**
     * @param _det the _det to set
     */
    public void setDet(int _det) {
        this._det = _det;
    }

    public Particle(int idx, int pid, double e, double px, double py, double pz, 
            double emc, double erec, double upx, double upy, double upz,
            double vx, double vy, double vz, 
            int charge,double mass,
            int[] ndau, int dau1idx, int dau2idx, int dau3idx) {
        this._Idx = idx;
        this._ndau = ndau;
        this._pid = pid;
        this._massConstrE = emc;
        this._E = e;
        this._uE = erec;
        this._px = px;
        this._py = py;
        this._pz = pz;
        this._vx = vx;
        this._vy = vy;
        this._vz = vz;
        this._upx = upx;
        this._upy = upy;
        this._upz = upz;
        this._charge = charge;
        this._mass = mass;
        //this._daughters = new ArrayList<Particle>();
        this.setP(Math.sqrt(_px*_px+_py*_py+_pz*_pz));
        this.setPt(Math.sqrt(_px*_px+_py*_py));
        this._daughters = new ArrayList<>(ndau.length);
//        for(int i =0; i< ndau; i++) {
//            this._daughters.add(new Particle());
//        }
//        if(ndau>1) {
//            this._daughters.get(0).setIdx(dau1idx);
//            this._daughters.get(1).setIdx(dau2idx);
//            if(ndau==3) {
//                this._daughters.get(2).setIdx(dau3idx);
//            } 
//        }
    }
    
    public boolean combine(Particle part1, Particle part2, int pid) {
        double mLo= Double.NEGATIVE_INFINITY;
        double mHi= Double.POSITIVE_INFINITY;
        
        return this.combine(part1, part2, pid, mLo, mHi);
    }
    
    public boolean combine(Particle part1, Particle part2, int pid, double mLo, double mHi) {
        double Vx = 0.5*(part1.getVx()+part2.getVx());
        double Vy = 0.5*(part1.getVy()+part2.getVy());
        double Vz = 0.5*(part1.getVz()+part2.getVz());
        
        double dx = part1.getVx()-part2.getVx();
        double dy = part1.getVy()-part2.getVy();
        double dz = part1.getVz()-part2.getVz();
        
        double r = Math.sqrt(dx*dx+dy*dy+dz*dz);
        
        double Px = part1.getPx()+part2.getPx();
        double Py = part1.getPy()+part2.getPy();
        double Pz = part1.getPz()+part2.getPz();
        
        double uPx = part1.getUpx()+part2.getUpx();
        double uPy = part1.getUpy()+part2.getUpy();
        double uPz = part1.getUpz()+part2.getUpz();
                
        double mass = this.calcMass(part1, part2);
        double uncormass = this.calcUncorrMass(part1, part2);
        double unconstmass = this.calcUnconstrainedMass(part1, part2);
        
        if(mass < mLo || mass > mHi) 
            return false;
        
        
        this.setPid(pid);
        this.setCharge(part1.getCharge()+part2.getCharge());
        this.setE(part1.getE()+part2.getE());
        this.setuE(part1.getuE()+part2.getuE());
        this.setMassConstrE(part1.getMassConstrE()+part2.getMassConstrE());
        this.setUmass(unconstmass);
        this.setUncormass(uncormass);
        this.setRecMass(mass);
        this.setMass(mass);
        this.setR(r);
        this.setVx(Vx);
        this.setVy(Vy);
        this.setVz(Vz);
        this.setPx(Px);
        this.setPy(Py);
        this.setPz(Pz);
        this.setUpx(uPx);
        this.setUpy(uPy);
        this.setUpz(uPz);
        this._daughters.add(part1);
        this._daughters.add(part2);
        //System.out.println("check"); System.out.println(part1.toString()); System.out.println(part2.toString());
        
        return true;
    }
    
    
    public double getECM() {
        double m = this.getMass();
        double px = this.getPxcm();
        double py = this.getPycm();
        double pz = this.getPzcm();
        
        return Math.sqrt(m*m+pz*px+py*py+pz*pz);
    }
    
    /**
     * @return the _Idx
     */
    public int getIdx() {
        return _Idx;
    }

    /**
     * @param _Idx the _Idx to set
     */
    public void setIdx(int _Idx) {
        this._Idx = _Idx;
    }
    /**
     * @return the _mass
     */
    public double getMass() {
        return _mass;
    }

    /**
     * @param _mass the _mass to set
     */
    public void setMass(double _mass) {
        this._mass = _mass;
    }

    /**
     * @return the _uncormass
     */
    public double getUncormass() {
        return _uncormass;
    }

    /**
     * @param _uncormass the _uncormass to set
     */
    public void setUncormass(double _uncormass) {
        this._uncormass = _uncormass;
    }
    
    /**
     * @return the _mass
     */
    public double getRecMass() {
        return _recmass;
    }

    /**
     * @param _recmass the _mass to set
     */
    public void setRecMass(double _mass) {
        this._recmass = _mass;
    }

    /**
     * @return the _pid
     */
    public int getPid() {
        return _pid;
    }

    /**
     * @param _pid the _pid to set
     */
    public void setPid(int _pid) {
        this._pid = _pid;
    }

    /**
     * @return the _E
     */
    public double getE() {
        return _E;
    }

    /**
     * @param _E the _E to set
     */
    public void setE(double _E) {
        this._E = _E;
    }

    /**
     * @return the _massConstrE
     */
    public double getMassConstrE() {
        return _massConstrE;
    }

    /**
     * @param _massConstrE the _massConstrE to set
     */
    public void setMassConstrE(double _massConstrE) {
        this._massConstrE = _massConstrE;
    }

    /**
     * @return the _px
     */
    public double getPx() {
        return _px;
    }

    /**
     * @param _px the _px to set
     */
    public void setPx(double _px) {
        this._px = _px;
    }

    /**
     * @return the _py
     */
    public double getPy() {
        return _py;
    }

    /**
     * @param _py the _py to set
     */
    public void setPy(double _py) {
        this._py = _py;
    }

    /**
     * @return the _pz
     */
    public double getPz() {
        return _pz;
    }

    /**
     * @param _pz the _pz to set
     */
    public void setPz(double _pz) {
        this._pz = _pz;
    }

    /**
     * @return the _vx
     */
    public double getVx() {
        return _vx;
    }

    /**
     * @param _vx the _vx to set
     */
    public void setVx(double _vx) {
        this._vx = _vx;
    }

    /**
     * @return the _vy
     */
    public double getVy() {
        return _vy;
    }

    /**
     * @param _vy the _vy to set
     */
    public void setVy(double _vy) {
        this._vy = _vy;
    }

    /**
     * @return the _vz
     */
    public double getVz() {
        return _vz;
    }

    /**
     * @param _vz the _vz to set
     */
    public void setVz(double _vz) {
        this._vz = _vz;
    }

    /**
     * @return the _r
     */
    public double getR() {
        return _r;
    }

    /**
     * @param _r the _r to set
     */
    public void setR(double _r) {
        this._r = _r;
    }

    /**
     * @return the _p
     */
    public double getP() {
        return _p;
    }

    /**
     * @param _p the _p to set
     */
    public void setP(double _p) {
        this._p = _p;
    }

    /**
     * @return the _pt
     */
    public double getPt() {
        return _pt;
    }

    /**
     * @param _pt the _pt to set
     */
    public void setPt(double _pt) {
        this._pt = _pt;
    }

    /**
     * @return the _pxcm
     */
    public double getPxcm() {
        return _pxcm;
    }

    /**
     * @param _pxcm the _pxcm to set
     */
    public void setPxcm(double _pxcm) {
        this._pxcm = _pxcm;
    }

    /**
     * @return the _pycm
     */
    public double getPycm() {
        return _pycm;
    }

    /**
     * @param _pycm the _pycm to set
     */
    public void setPycm(double _pycm) {
        this._pycm = _pycm;
    }

    /**
     * @return the _pzcm
     */
    public double getPzcm() {
        return _pzcm;
    }

    /**
     * @param _pzcm the _pzcm to set
     */
    public void setPzcm(double _pzcm) {
        this._pzcm = _pzcm;
    }

    /**
     * @return the _beta
     */
    public double getBeta() {
        return _beta;
    }

    /**
     * @param _beta the _beta to set
     */
    public void setBeta(double _beta) {
        this._beta = _beta;
    }

    /**
     * @return the _calcBeta
     */
    public double getCalcBeta() {
        return _calcBeta;
    }

    /**
     * @param _calcBeta the _calcBeta to set
     */
    public void setCalcBeta(double _calcBeta) {
        this._calcBeta = _calcBeta;
    }

    /**
     * @return the _charge
     */
    public int getCharge() {
        return _charge;
    }

    /**
     * @param _charge the _charge to set
     */
    public void setCharge(int _charge) {
        this._charge = _charge;
    }

    
    /**
     * @return the _status
     */
    public int getStatus() {
        return _status;
    }

    /**
     * @param _status the _status to set
     */
    public void setStatus(int _status) {
        this._status = _status;
    }

    /**
     * @return the _daughters
     */
    public List<Particle> getDaughters() {
        return _daughters;
    }

    /**
     * @param _daughters the _daughters to set
     */
    public void setDaughters(List<Particle> _daughters) {
        this._daughters = _daughters;
    }
    
    public boolean isFT() {
        return Math.abs(this.getStatus())/1000==1;
    }

    public boolean isFD() {
        return Math.abs(this.getStatus())/1000==2;
    }

    public boolean isCD() {
        return Math.abs(this.getStatus())/1000==4;
    }

    private double calcMass(Particle part1, Particle part2) {
        double E1 = part1.getMassConstrE();
        double E2 = part2.getMassConstrE();
        //double E1 = part1.getE();
        //double E2 = part2.getE();
        double Px = part1.getPx()+part2.getPx();
        double Py = part1.getPy()+part2.getPy();
        double Pz = part1.getPz()+part2.getPz();
        
        double E = E1+E2;
        return Math.sqrt(E*E - Px*Px-Py*Py-Pz*Pz);
        
    }
    private double calcUnconstrainedMass(Particle part1, Particle part2) {
        double E1 = part1.getE();
        double E2 = part2.getE();
        double Px = part1.getPx()+part2.getPx();
        double Py = part1.getPy()+part2.getPy();
        double Pz = part1.getPz()+part2.getPz();
        
        double E = E1+E2;
        return Math.sqrt(E*E - Px*Px-Py*Py-Pz*Pz);
        
    }
    private double calcUncorrMass(Particle part1, Particle part2) { //momentum not corrected at the vertex
        double E1 = part1.getMassConstrE();
        double E2 = part2.getMassConstrE();
        //double E1 = part1.getE();
        //double E2 = part2.getE();
        double Px = part1.getUpx()+part2.getUpx();
        double Py = part1.getUpy()+part2.getUpy();
        double Pz = part1.getUpz()+part2.getUpz();
        
        double E = E1+E2;
        return Math.sqrt(E*E - Px*Px-Py*Py-Pz*Pz);
        
    }
    private Swim swim;
    double[][] SwimTo(Particle p2, Swim swim) {
        double vx1,vy1,vz1,vx2,vy2,vz2;
        double px1,py1,pz1,px2,py2,pz2;
        int q = this.getCharge(); 
        if(this.getVx()==0 && this.getVy()==0 && this.getVz()==0) {
            vx1 = this.getOvx();
            vy1 = this.getOvy();
            vz1 = this.getOvz();
            px1 = this.getUpx();
            py1 = this.getUpy();
            pz1 = this.getUpz();
        } else {
            vx1 = this.getVx();
            vy1 = this.getVy();
            vz1 = this.getVz();
            px1 = this.getPx();
            py1 = this.getPy();
            pz1 = this.getPz();
        }
        if(p2.getVx()==0 && p2.getVy()==0 && p2.getVz()==0) {
            vx2 = p2.getOvx();
            vy2 = p2.getOvy();
            vz2 = p2.getOvz();
            px2 = p2.getUpx();
            py2 = p2.getUpy();
            pz2 = p2.getUpz();
        } else {
            vx2 = p2.getVx();
            vy2 = p2.getVy();
            vz2 = p2.getVz();
            px2 = p2.getPx();
            py2 = p2.getPy();
            pz2 = p2.getPz();
        }
        
        return this.SwimTo(q, vx1, vy1, vz1, px1, py1, pz1, vx2, vy2, vz2, px2, py2, pz2, swim);
    }  
    double[][] SwimTo(int q, double vxch, double vych, double vzch, double pxch, double pych, double pzch,
            double vxvo, double vyvo, double vzvo, double pxvo, double pyvo, double pzvo, Swim swim) {
        
        double[][] result = new double[2][6]; //[0] the charged track; [1] the straight track
        if(vzvo>vzch) {
            swim.SetSwimParameters(vxch, vych, vzch, 
                pxch, pych, pzch, q);
        } else {
            swim.SetSwimParameters(vxch, vych, vzch, 
                -pxch, -pych, -pzch, -q);
        }
        double p = Math.sqrt(pxvo*pxvo+pyvo*pyvo+pzvo*pzvo);
        double theta = Math.toDegrees(Math.acos(pzvo/p));
        double phi = Math.toDegrees(Math.atan2(pyvo, pxvo));
        
        double step = 0.02; //cm
        double vxt = vxvo/100.0; //swimTrajctory arg in m
        double vyt = vyvo/100.0;
        double vzt = vzvo/100.0;
        
        
        int nstep = (int)(Math.sqrt((vxvo-vxch)*(vxvo-vxch)+(vyvo-vych)*(vyvo-vych)+(vzvo-vzch)*(vzvo-vzch))/step)+2;
        
        SwimTrajectory tr = new SwimTrajectory(q,vxt-nstep*pxvo/p*step/100.0,vyt-nstep*pyvo/p*step/100.0,vzt-nstep*pzvo/p*step/100.0,p,theta,phi,nstep+2);
        for(int i =nstep-1; i>-1; i--) {
            tr.add(vxt-i*pxvo/p*step/100.0,vyt-i*pyvo/p*step/100.0,vzt-i*pzvo/p*step/100.0, p, theta, phi); 
        }
        double[] pars = swim.SwimToDCA(tr);
        result[0][0] = pars[0];
        result[0][1] = pars[1];
        result[0][2] = pars[2];
        result[0][3] = pars[3];
        result[0][4] = pars[4];
        result[0][5] = pars[5];
        
        Line3D l = new Line3D(vxvo-100*pxvo/p, vyvo-100*pyvo/p, vzvo-100*pzvo/p,
                                vxvo+100*pxvo/p, vyvo+100*pyvo/p, vzvo+100*pzvo/p);
        
        Point3D p2 = l.distance(new Point3D(pars[0], pars[1], pars[2])).origin();
        result[1][0] = p2.x();
        result[1][1] = p2.y();
        result[1][2] = p2.z();
        result[1][3] = pxvo;
        result[1][4] = pyvo;
        result[1][5] = pzvo;
        
        return result;
    }
    public boolean hasDaughters() {
        boolean d = false;
        if(this.getDaughters()!=null && !this.getDaughters().isEmpty())
            d = true;
        return d;
    }
    private String printVP(Particle p) {
        String s = "\n Dec X ";
        s+=new Point3D(p.getVx(), p.getVy(), p.getVz()).toString()+ "\n Prod X "+
                new Point3D(p.getOvx(), p.getOvy(), p.getOvz()).toString()
                +"\n P "+new Point3D(p.getPx(), p.getPy(), p.getPz()).toString()
                +"\n uP "+new Point3D(p.getUpx(), p.getUpy(), p.getUpz()).toString();
        return s;
    }
    @Override
    public String toString() {
        String s="";
        s+=this.getIdx()+" pid "+this.getPid()+printVP(this)+"\n mass "+this.getMass()+"\n ndau "+this.getDaughters().size();
        if(this.hasDaughters()) {
            for(Particle p : this.getDaughters()) {
                if(p!=null) {
                    s+="\ndaughter "+p.getIdx()+" pid "+p.getPid()+printVP(p)+"\n mass "+p.getMass();
                    if(p.hasDaughters()) {
                        for(int i = 0; i<p.getDaughters().size(); i++) {
                            Particle p2 = p.getDaughters().get(i);
                            s+="\n granddaughter "+p2.getIdx()+" pid "+p2.getPid()+printVP(p2)+"\n mass "+p2.getMass();
                        }
                    }
                }
            }
        } 
        return s;
    }

    
}
