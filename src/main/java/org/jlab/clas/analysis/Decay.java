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
import org.jlab.io.base.DataBank;

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
            List<Particle> daughters) {
        _parPID = parPID;
        _dau1PID = dau1PID;
        _dau2PID = dau2PID;
        _dau3PID = dau3PID;
        _loMassCut = loMassCut;
        _hiMassCut = hiMassCut;
        _daughters = daughters;
        
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
                            if(part1.getCharge()!=0 && part2.getCharge()!=0) {
                                part1.setVx(vx1);
                                part1.setVy(vy1);
                                part1.setVz(vz1);
                                part1.setPx(px1);
                                part1.setPy(py1);
                                part1.setPz(pz1);
                                part2.setVx(vx2);
                                part2.setVy(vy2);
                                part2.setVz(vz2);
                                part2.setPx(px2);
                                part2.setPy(py2);
                                part2.setPz(pz2);
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
    
    private double r = 999;
    private double vx =999;
    private double vy =999;
    private double vz =999;
    private double vx1 =999;
    private double vy1 =999;
    private double vz1 =999;
    private double px1 =999;
    private double py1 =999;
    private double pz1 =999;
    private double vx2 =999;
    private double vy2 =999;
    private double vz2 =999;
    private double px2 =999;
    private double py2 =999;
    private double pz2 =999;
    
    private void reset() {
        r = 999;
        vx =999;
        vy =999;
        vz =999;
        vx1 =999;
        vy1 =999;
        vz1 =999;
        px1 =999;
        py1 =999;
        pz1 =999;
        vx2 =999;
        vy2 =999;
        vz2 =999;
        px2 =999;
        py2 =999;
        pz2 =999;
    }
    private boolean checkVertex(Particle p1, Particle p2) {
        if(getVertBank()==null) return false;
        if(getVertBank()!=null) {
            int nrows2 = getVertBank().rows();
            for(int loop2 = 0; loop2 < nrows2; loop2++){
                reset();
                int index1 = -1;//(int) getVertBank().getShort("index1", loop2);
                int index2 = -1;//(int) getVertBank().getShort("index2", loop2);
                if(p1.getIdx()-1==(int) getVertBank().getShort("index1", loop2)) {
                    index1 = (int) getVertBank().getShort("index1", loop2);
                    if(p2.getIdx()-1==(int) getVertBank().getShort("index2", loop2)) {
                        index2 = (int) getVertBank().getShort("index2", loop2);
                    }
                    if(index1!=-1 && index2!=-1) {
                        p1.vIndex=loop2;
                        p2.vIndex=loop2;
                        r =  (double) getVertBank().getFloat("r", loop2);
                        vx = (double) getVertBank().getFloat("x", loop2);
                        vy = (double) getVertBank().getFloat("y", loop2);
                        vz = (double) getVertBank().getFloat("z", loop2);
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
                    }
                } 
                if(p1.getIdx()-1==(int) getVertBank().getShort("index2", loop2)) {
                    index1 = (int) getVertBank().getShort("index2", loop2);
                    if(p2.getIdx()-1==(int) getVertBank().getShort("index1", loop2)) {
                        index2 = (int) getVertBank().getShort("index1", loop2);
                    }
                    if(index1!=-1 && index2!=-1) {
                        p1.vIndex=loop2;
                        p2.vIndex=loop2;
                        r =  (double) getVertBank().getFloat("r", loop2);
                        vx = (double) getVertBank().getFloat("x", loop2);
                        vy = (double) getVertBank().getFloat("y", loop2);
                        vz = (double) getVertBank().getFloat("z", loop2);
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
                    }
                }
            }
            if(r!=999) {
                return true;
            } 
        }
        return false;
    }
}
