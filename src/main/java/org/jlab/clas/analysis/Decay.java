/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.swimtools.Swim;
import org.jlab.io.base.DataBank;
import org.jlab.rec.vtx.VertexFinder;

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
    private boolean useSwimmer = true;
    
    public Decay(int parPID, int dau1PID, int dau2PID, int dau3PID, double loMassCut, double hiMassCut,
            List<Particle> daughters, double xb, double yb, Swim swimmer) {
        this.setSwimmer(swimmer);
        _parPID = parPID;
        _dau1PID = dau1PID;
        _dau2PID = dau2PID;
        _dau3PID = dau3PID;
        _loMassCut = loMassCut;
        _hiMassCut = hiMassCut;
        _daughters = daughters;
        
        this.setxB(xb);
        this.setyB(yb);
        //System.out.println("PIDS "+dau1PID+", "+dau2PID);
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
                    for(Particle part1 : list1) {
                        for(Particle part2 : list2) { 
                            Particle part = new Particle();
                            if(part1.getCharge()!=0 && part2.getCharge()!=0 ) {
                                    if(this.checkVertex(part1,part2)==false)
                                        continue;
                            }
                            if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)) { 
                                
                                if(part1.getCharge()!=0 && part2.getCharge()==0) {
                                    vx = part1.getVx();
                                    vy = part1.getVy();
                                    vz = part1.getVz();
                                } else if(part2.getCharge()!=0 && part1.getCharge()==0) {
                                    vx = part2.getVx();
                                    vy = part2.getVy();
                                    vz = part2.getVz();
                                } else if(part2.getCharge()!=0 && part1.getCharge()==0) {
                                    vx = 0;
                                    vy = 0;
                                    vz = 0;
                                }
                                part1.isUsed = true;
                                part2.isUsed = true;
                                part.setVx(vx);
                                part.setVy(vy);
                                part.setVz(vz);
                                part.setR(r);
                                _particles.add(part);
                            }
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

    private void flagCombinatorials() {
        if(r==999) 
            return;

        int j = 99;
        for(Particle p : _particles) {
            p.setIdx(j++);
        }
        for(Particle p1 : _particles) {
            for(Particle p2 : _particles) {
                if(p1.getIdx()==p2.getIdx())
                    continue;
                if(p1.getDaughters().get(0).getIdx()==p2.getDaughters().get(0).getIdx() || 
                        p1.getDaughters().get(1).getIdx()==p2.getDaughters().get(1).getIdx() ) {
                    this.reject(p1,p2).setIdx(-1);
                }
            }
        }
       
    }

    private Particle reject(Particle p1, Particle p2) {
        double mass = PDGDatabase.getParticleById(Math.abs(p1.getPid())).mass();
        if(Math.abs(p1.getRecMass()-mass)>Math.abs(p2.getRecMass()-mass)) {
            return p1;
        } else {
            return p2;
        }
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
    private Swim swimmer;
    private double xB;
    private double yB;
    
    private double r = 999;
    private double vx =999;
    private double vy =999;
    private double vz =999;
    
    private boolean checkVertex(Particle p1, Particle p2) {
        if(getVertBank()!=null) {
            int nrows2 = getVertBank().rows();
            for(int loop2 = 0; loop2 < nrows2; loop2++){
                int index1 = (int) getVertBank().getShort("index1", loop2);
                int index2 = (int) getVertBank().getShort("index2", loop2);

                if(p1.getIdx()-1==index1 || p1.getIdx()-1==index2)
                    p1.vIndex=loop2;
                if(p2.getIdx()-1==index1 || p2.getIdx()-1==index2)
                    p2.vIndex=loop2;
                if(p1.vIndex==p2.vIndex) 
                    loop2 = nrows2;
            }
            if(p1.vIndex==-1 || p2.vIndex==-1 || p1.vIndex!=p2.vIndex) 
                return false;
            r =  (double) getVertBank().getFloat("r", p1.vIndex);
            vx = (double) getVertBank().getFloat("x", p1.vIndex);
            vy = (double) getVertBank().getFloat("y", p1.vIndex);
            vz = (double) getVertBank().getFloat("z", p1.vIndex);
            
            return true;
        } else {
            return false;
        }
    }

    /**
     * @return the swimmer
     */
    public Swim getSwimmer() {
        return swimmer;
    }

    /**
     * @param swimmer the swimmer to set
     */
    public void setSwimmer(Swim swimmer) {
        this.swimmer = swimmer;
    }

    /**
     * @return the xB
     */
    public double getxB() {
        return xB;
    }

    /**
     * @param xB the xB to set
     */
    public void setxB(double xB) {
        this.xB = xB;
    }

    /**
     * @return the yB
     */
    public double getyB() {
        return yB;
    }

    /**
     * @param yB the yB to set
     */
    public void setyB(double yB) {
        this.yB = yB;
    }
    
    VertexFinder vertexFinder = new VertexFinder();

        
    
}
