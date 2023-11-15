/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.clas.analysis.event.Reader;
import org.jlab.clas.swimtools.Swim;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.io.base.DataBank;
import org.jlab.rec.vtx.DoubleSwim;

/**
 *
 * @author ziegler
 */
public class Decay extends Particle {
    
    private int _parPID;
    private int _dau1PID;
    private int _dau2PID;
    private int _dau3PID;
    
    private double _loMassCut;
    private double _hiMassCut;
    
    private List<Particle> _daughters;
    private List<Particle> _particles;
    
    private Map<Integer, List<Particle>> listsByPID = new HashMap<>();
    private static DataBank vertBank;
    
    public Decay(int parPID, int dau1PID, int dau2PID, int dau3PID, double loMassCut, double hiMassCut, 
            double loMMassCut, double hiMMassCut, double eBeam, Reader.Target target,  
            List<Particle> daughters, Particle el, Swim swim, int pass) { 
        _parPID = parPID;
        _dau1PID = dau1PID;
        _dau2PID = dau2PID;
        _dau3PID = dau3PID;
        _loMassCut = loMassCut;
        _hiMassCut = hiMassCut;
        _daughters = daughters;
        
        this.sortByPID(daughters, dau1PID, dau2PID, dau3PID);
        List<Particle> list1 = new ArrayList<>();
        List<Particle> list2 = new ArrayList<>();
        List<Particle> list3 = new ArrayList<>();
        
        if(listsByPID!=null) { 
            if(listsByPID.containsKey(dau1PID) && dau1PID!=0)
                list1.addAll(listsByPID.get(dau1PID));
            if(listsByPID.containsKey(dau2PID) && dau2PID!=0)
                list2.addAll(listsByPID.get(dau2PID));
            if(listsByPID.containsKey(dau3PID) && dau3PID!=0)
                list3.addAll(listsByPID.get(dau3PID)); 
            if(!list1.isEmpty() && !list2.isEmpty()) {
                _particles = new ArrayList<>();
                if(list3.isEmpty()) {
                    for(Particle part01 : list1) {
                        for(Particle part02 : list2) { 
                            Particle part = new Particle();
                            boolean overlaps = this.overlaps(part01, part02);
                            if(overlaps) { 
                                continue;
                            }
                            //make a clone and assign vtx row id
                            Particle part1 = new Particle(part01);
                            Particle part2 = new Particle(part02); 
                            if(part1.getCharge()!=0 && part2.getCharge()!=0 ) { 
                                if(pass==1)
                                    if(this.checkVertex(part1,part2)==false) { 
                                       continue;
                                    } 
                                if(pass>1 )
                                    if(this.makeVertex(part1,part2)==false) {
                                        continue;
                                    }
                                if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)==false) 
                                    continue;
                            }
                            
                            if(part1.getCharge()!=0 && part2.getCharge()==0) { 
                                
                                double[] p1 = part1.SwimTo(part2, swim);
                                double v1x = p1[0];
                                double v1y = p1[1];
                                double v1z = p1[2];
                                double p1x = p1[3];
                                double p1y = p1[4];
                                double p1z = p1[5];
                                int q1 = part1.getCharge();
                                double ue1 = part1.getuE();
                                double e1 = part1.getE();
                                double emc1 = part1.getMassConstrE();
                                double up1x = part1.getUpx();
                                double up1y = part1.getUpy();
                                double up1z = part1.getUpz();
                                int q2 = part2.getCharge();
                                double ue2 = part2.getuE();
                                double e2 = part2.getE();
                                double emc2 = part2.getMassConstrE();
                                double up2x = part2.getUpx();
                                double up2y = part2.getUpy();
                                double up2z = part2.getUpz();
                                
                                double v2x = part2.getVx();
                                double v2y = part2.getVy();
                                double v2z = part2.getVz();
                                double p2x = part2.getPx();
                                double p2y = part2.getPy();
                                double p2z = part2.getPz();
                                
                                if(part.combine(q1, ue1, e1, emc1, v1x, v1y, v1z, p1x, p1y, p1z, up1x, up1y, up1z, 
                                        q2, ue2, e2, emc2, v2x, v2y, v2z, p2x, p2y, p2z, up2x, up2y, up2z, parPID, loMassCut, hiMassCut)==false) continue;
                                part1.setVx(v1x);
                                part1.setVy(v1y);
                                part1.setVz(v1z);
                                part1.setPx(p1x);
                                part1.setPy(p1y);
                                part1.setPz(p1z);
                                part.getDaughters().add(part1);
                                part.getDaughters().add(part2);
                                
                            } else if(part2.getCharge()!=0 && part1.getCharge()==0) { 
                                double[] p2 = part2.SwimTo(part1, swim);
                                double v2x = p2[0];
                                double v2y = p2[1];
                                double v2z = p2[2];
                                double p2x = p2[3];
                                double p2y = p2[4];
                                double p2z = p2[5];
                                int q1 = part1.getCharge();
                                double ue1 = part1.getuE();
                                double e1 = part1.getE();
                                double emc1 = part1.getMassConstrE();
                                double up1x = part1.getUpx();
                                double up1y = part1.getUpy();
                                double up1z = part1.getUpz();
                                int q2 = part2.getCharge();
                                double ue2 = part2.getuE();
                                double e2 = part2.getE();
                                double emc2 = part2.getMassConstrE();
                                double up2x = part2.getUpx();
                                double up2y = part2.getUpy();
                                double up2z = part2.getUpz();
                                
                                double v1x = part1.getVx();
                                double v1y = part1.getVy();
                                double v1z = part1.getVz();
                                double p1x = part1.getPx();
                                double p1y = part1.getPy();
                                double p1z = part1.getPz();
                                
                                if(part.combine(q1, ue1, e1, emc1, v1x, v1y, v1z, p1x, p1y, p1z, up1x, up1y, up1z, 
                                        q2, ue2, e2, emc2, v2x, v2y, v2z, p2x, p2y, p2z, up2x, up2y, up2z, parPID, loMassCut, hiMassCut)==false) continue;
                                part2.setVx(v2x);
                                part2.setVy(v2y);
                                part2.setVz(v2z);
                                part2.setPx(p2x);
                                part2.setPy(p2y);
                                part2.setPz(p2z);
                                part.getDaughters().add(part1);
                                part.getDaughters().add(part2);
                            } else if(part2.getCharge()==0 && part1.getCharge()==0) { 
                                Point3D P1X1 = new Point3D(part1.getVx()-100*part1.getPx(), 
                                        part1.getVy()-100*part1.getPy(), 
                                        part1.getVz()-100*part1.getPz());
                                Point3D P1X2 = new Point3D(part1.getVx()+100*part1.getPx(), 
                                        part1.getVy()+100*part1.getPy(), 
                                        part1.getVz()+100*part1.getPz());
                                Point3D P2X1 = new Point3D(part2.getVx()-100*part2.getPx(), 
                                        part2.getVy()-100*part2.getPy(), 
                                        part2.getVz()-100*part2.getPz());
                                Point3D P2X2 = new Point3D(part2.getVx()+100*part2.getPx(), 
                                        part2.getVy()+100*part2.getPy(), 
                                        part2.getVz()+100*part2.getPz());
                                Line3D par1Extrap =new Line3D(P1X1,P1X2);
                                Line3D par2Extrap =new Line3D(P2X1,P2X2);
                                Line3D parExtrap = par1Extrap.distance(par2Extrap);
                                double v1x = parExtrap.origin().x();
                                double v1y = parExtrap.origin().y();
                                double v1z = parExtrap.origin().z();
                                double v2x = parExtrap.end().x();
                                double v2y = parExtrap.end().y();
                                double v2z = parExtrap.end().z();
                                
                                int q1 = part1.getCharge();
                                double ue1 = part1.getuE();
                                double e1 = part1.getE();
                                double emc1 = part1.getMassConstrE();
                                double p1x = part1.getPx();
                                double p1y = part1.getPy();
                                double p1z = part1.getPz();
                                double up1x = part1.getUpx();
                                double up1y = part1.getUpy();
                                double up1z = part1.getUpz();
                                int q2 = part2.getCharge();
                                double ue2 = part2.getuE();
                                double e2 = part2.getE();
                                double emc2 = part2.getMassConstrE();
                                double p2x = part2.getPx();
                                double p2y = part2.getPy();
                                double p2z = part2.getPz();
                                double up2x = part2.getUpx();
                                double up2y = part2.getUpy();
                                double up2z = part2.getUpz();
                               
                                if(part.combine(q1, ue1, e1, emc1, v1x, v1y, v1z, p1x, p1y, p1z, up1x, up1y, up1z, 
                                                q2, ue2, e2, emc2, v2x, v2y, v2z, p2x, p2y, p2z, up2x, up2y, up2z, 
                                                parPID, loMassCut, hiMassCut)==false) continue;
                                part.getDaughters().add(part1);
                                part.getDaughters().add(part2);
                               
                            }
                            part1.isUsed = true;
                            part2.isUsed = true;
                            
                            if(part.getMissingMass(el, eBeam, target)>loMMassCut && part.getMissingMass(el, eBeam, target)<hiMMassCut)
                                _particles.add(part);
                        
                        }
                    }
                    //this.flagCombinatorials();
                } else {
                    for(Particle part1 : list1) {
                        for(Particle part2 : list2) {
                            for(Particle part3 : list3) {
                                Particle part = new Particle();
                            //if(part.combine(part1, part2, part3, parPID, hiMassCut, hiMassCut))
                            //    _particles.add(part);
                            }
                        }
                    }
                }
            }
        }
    }

    private void sortByPID(List<Particle> daughters, int dau1PID, int dau2PID, int dau3PID) {
        listsByPID.clear();
        for(Particle par : daughters) { 
            if(par.getPid()==dau1PID || par.getPid()==dau2PID || par.getPid()==dau3PID) {
                if(listsByPID.containsKey(par.getPid())) {
                    listsByPID.get(par.getPid()).add(par); 
                } else {
                    listsByPID.put(par.getPid(), new ArrayList<Particle>());
                    listsByPID.get(par.getPid()).add(par);  
                }
            }
        }
    }

    /**
     * @return the _particles
     */
    public List<Particle> getParticles() {
        return _particles;
    }

    /**
     * @param _particles the _particles to set
     */
    public void setParticles(List<Particle> _particles) {
        this._particles = _particles;
    }

    /**
     * @return the vertBank
     */
    public DataBank getVertBank() {
        return vertBank;
    }

    /**
     * @param vertBank the vertBank to set
     */
    public static void setVertBank(DataBank vertBank) {
        Decay.vertBank = vertBank;
    }
    private DoubleSwim _ds = new DoubleSwim();
    private boolean makeVertex(Particle p1, Particle p2) { 
        boolean pass = false;
        double vx1,vy1,vz1,vx2,vy2,vz2;
        double px1,py1,pz1,px2,py2,pz2; 
        if(p1.getVx()==0 && p1.getVy()==0 && p1.getVz()==0) {
            vx1 = p1.getOvx();
            vy1 = p1.getOvy();
            vz1 = p1.getOvz();
            px1 = p1.getUpx();
            py1 = p1.getUpy();
            pz1 = p1.getUpz();
        } else {
            vx1 = p1.getVx();
            vy1 = p1.getVy();
            vz1 = p1.getVz();
            px1 = p1.getPx();
            py1 = p1.getPy();
            pz1 = p1.getPz();
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
        _ds.init(vx1,vy1,vz1,px1,py1,pz1, p1.getCharge(),
                 vx2,vy2,vz2,px2,py2,pz2, p2.getCharge());

        double[][] t = _ds.getDoubleSwimVertexes();
        
        if(t==null) return pass;
            
        double r = Math.sqrt((t[0][0]-t[1][0])*(t[0][0]-t[1][0])
                     +(t[0][1]-t[1][1])*(t[0][1]-t[1][1])
                     +(t[0][2]-t[1][2])*(t[0][2]-t[1][2])  );

        p1.setVx(t[0][0]);
        p1.setVy(t[0][1]);
        p1.setVz(t[0][2]);
        p1.setPx(t[0][3]);
        p1.setPy(t[0][4]);
        p1.setPz(t[0][5]);
        p2.setVx(t[1][0]);
        p2.setVy(t[1][1]);
        p2.setVz(t[1][2]);
        p2.setPx(t[1][3]);
        p2.setPy(t[1][4]);
        p2.setPz(t[1][5]);
        
        pass= true;
        
        return pass;
     }
    
    private boolean checkVertDocaBank(Particle p1, Particle p2) {
        boolean pass = false;
        if(getVertBank()==null) return pass;
        
        if(getVertBank()!=null) { 
            int nrows2 = getVertBank().rows();
            for(int loop2 = 0; loop2 < nrows2; loop2++){
                double r = 999;
                //double vx =999;
                //double vy =999;
                //double vz =999;
                double vx1 =999;
                double vy1 =999;
                double vz1 =999;
                double px1 =999;
                double py1 =999;
                double pz1 =999;
                double vx2 =999;
                double vy2 =999;
                double vz2 =999;
                double px2 =999;
                double py2 =999;
                double pz2 =999;
                
                if(p1.getIdx()-1==(int) getVertBank().getShort("index1", loop2)
                        && p2.getIdx()-1==(int) getVertBank().getShort("index2", loop2)) {
                    p1.vIndex=loop2;
                    p2.vIndex=loop2;
                    
                    r =  (double) getVertBank().getFloat("r", loop2);
                    //vx = (double) getVertBank().getFloat("x", loop2);
                    //vy = (double) getVertBank().getFloat("y", loop2);
                    //vz = (double) getVertBank().getFloat("z", loop2);
                    vx1 = (double) getVertBank().getFloat("x1", loop2);
                    vy1 = (double) getVertBank().getFloat("y1", loop2);
                    vz1 = (double) getVertBank().getFloat("z1", loop2);
                    px1 = (double) getVertBank().getFloat("cx1", loop2);
                    py1 = (double) getVertBank().getFloat("cy1", loop2);
                    pz1 = (double) getVertBank().getFloat("cz1", loop2);
                    vx2 = (double) getVertBank().getFloat("x2", loop2);
                    vy2 = (double) getVertBank().getFloat("y2", loop2);
                    vz2 = (double) getVertBank().getFloat("z2", loop2);
                    px2 = (double) getVertBank().getFloat("cx2", loop2);
                    py2 = (double) getVertBank().getFloat("cy2", loop2);
                    pz2 = (double) getVertBank().getFloat("cz2", loop2);
                    p1.setVx(vx1);
                    p1.setVy(vy1);
                    p1.setVz(vz1);
                    p1.setPx(px1);
                    p1.setPy(py1);
                    p1.setPz(pz1);
                    p2.setVx(vx2);
                    p2.setVy(vy2);
                    p2.setVz(vz2);
                    p2.setPx(px2);
                    p2.setPy(py2);
                    p2.setPz(pz2);
                    pass=true;
                    return pass;
                }
                if(p2.getIdx()-1==(int) getVertBank().getShort("index1", loop2)
                        && p1.getIdx()-1==(int) getVertBank().getShort("index2", loop2)) {
                    p1.vIndex=loop2;
                    p2.vIndex=loop2;
                    //r =  (double) getVertBank().getFloat("r", loop2);
                    //vx = (double) getVertBank().getFloat("x", loop2);
                    //vy = (double) getVertBank().getFloat("y", loop2);
                    //vz = (double) getVertBank().getFloat("z", loop2);
                    vx1 = (double) getVertBank().getFloat("x2", loop2);
                    vy1 = (double) getVertBank().getFloat("y2", loop2);
                    vz1 = (double) getVertBank().getFloat("z2", loop2);
                    px1 = (double) getVertBank().getFloat("cx2", loop2);
                    py1 = (double) getVertBank().getFloat("cy2", loop2);
                    pz1 = (double) getVertBank().getFloat("cz2", loop2);
                    vx2 = (double) getVertBank().getFloat("x1", loop2);
                    vy2 = (double) getVertBank().getFloat("y1", loop2);
                    vz2 = (double) getVertBank().getFloat("z1", loop2);
                    px2 = (double) getVertBank().getFloat("cx1", loop2);
                    py2 = (double) getVertBank().getFloat("cy1", loop2);
                    pz2 = (double) getVertBank().getFloat("cz1", loop2);
                    p1.setVx(vx1);
                    p1.setVy(vy1);
                    p1.setVz(vz1);
                    p1.setPx(px1);
                    p1.setPy(py1);
                    p1.setPz(pz1);
                    p2.setVx(vx2);
                    p2.setVy(vy2);
                    p2.setVz(vz2);
                    p2.setPx(px2);
                    p2.setPy(py2);
                    p2.setPz(pz2);
                    //this.resetE(p1,p2);
                    pass=true; 
                    return pass;
                }
                
            }
            
        }
       
        return pass;
    }    
    
    
    private boolean checkVertex(Particle p1, Particle p2) {
        boolean pass = this.checkVertDocaBank(p1, p2);
        if(getVertBank()==null) 
            pass = this.makeVertex(p1, p2);
        
        return pass;
    }

    private boolean overlaps(Particle part01, Particle part02) {
        boolean overlaps = false;
        
        List<Particle> part1 = new ArrayList<>();
        List<Particle> part2 = new ArrayList<>();
        part1.add(part01);
        part2.add(part02);
        
        for(int i1 = 0; i1<part01.getDaughters().size(); i1++) {
            part1.add(part01.getDaughters().get(i1)); 
            for(int ii1 = 0; ii1<part01.getDaughters().get(i1).getDaughters().size(); ii1++) {
                part1.add(part01.getDaughters().get(i1).getDaughters().get(ii1));
                for(int iii1 = 0; iii1<part01.getDaughters().get(i1).getDaughters().get(ii1).getDaughters().size(); iii1++) {
                    part1.add(part01.getDaughters().get(i1).getDaughters().get(ii1).getDaughters().get(iii1));
                }
            }
        }
        for(int i2 = 0; i2<part02.getDaughters().size(); i2++) {
            part2.add(part02.getDaughters().get(i2));
            for(int ii2 = 0; ii2<part02.getDaughters().get(i2).getDaughters().size(); ii2++) {
                part2.add(part02.getDaughters().get(i2).getDaughters().get(ii2));
                for(int iii2 = 0; iii2<part02.getDaughters().get(i2).getDaughters().get(ii2).getDaughters().size(); iii2++) {
                    part2.add(part02.getDaughters().get(i2).getDaughters().get(ii2).getDaughters().get(iii2));
                }
            }
        }
        for(Particle p1 : part1) {
            for(Particle p2 : part2) { 
                if(p1.getIdx()==p2.getIdx()) {
                    overlaps = true;
                }
            }
        }
        
        return overlaps;
    }

    private void resetE(Particle p1, Particle p2) {
        double beta1 = p1.getBeta();
        double p_1 = Math.sqrt(p1.getPx()*p1.getPx()+p1.getPy()*p1.getPy()+p1.getPz()*p1.getPz());
        double m_1 = p1.getMass();
        p1.setUncormass(Math.sqrt((p_1/beta1)*(p_1/beta1)-p_1*p_1));
        double mu_1 = p1.getUncormass();
        double beta2 = p2.getBeta();
        double p_2 = Math.sqrt(p2.getPx()*p2.getPx()+p2.getPy()*p2.getPy()+p2.getPz()*p2.getPz());
        p2.setUncormass(Math.sqrt((p_2/beta2)*(p_2/beta2)-p_2*p_2));
        double m_2 = p2.getMass();
        double mu_2 = p2.getUncormass();
        double E1 = Math.sqrt(p_1*p_1+m_1*m_1);
        double E2 = Math.sqrt(p_2*p_2+m_2*m_2);
        double uE1 = Math.sqrt(p_1*p_1+mu_1*mu_1);
        double uE2 = Math.sqrt(p_2*p_2+mu_2*mu_2);
        p1.setE(E1);
        p2.setE(E2);
        p1.setMassConstrE(E1);
        p2.setMassConstrE(E2);
        p1.setuE(uE1);
        p2.setuE(uE2); 
    }
    
}
