/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
            List<Particle> daughters, Swim swim, int pass) { 
        System.out.println("REACTION "+parPID+"-->"+dau1PID+" "+dau2PID);
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
                        /*    if(part01.getDaughters().size()==2 && part02.getDaughters().size()==2) {
                                if(part01.getDaughters().get(0).getIdx()==part02.getDaughters().get(0).getIdx() || 
                                        part01.getDaughters().get(1).getIdx()==part02.getDaughters().get(0).getIdx() || 
                                        part01.getDaughters().get(0).getIdx()==part02.getDaughters().get(1).getIdx() || 
                                        part01.getDaughters().get(1).getIdx()==part02.getDaughters().get(1).getIdx() )
                                    continue;
                            }
                         */
                            System.out.println(part01.toString());
                            System.out.println(part02.toString());
                            if(overlaps) { 
                                continue;
                            }
                            System.out.println("no overlaps "+part02.getCharge());
                            //make a clone and assign vtx row id
                            Particle part1 = new Particle(part01);
                            Particle part2 = new Particle(part02);
                            if(part1.getCharge()!=0 && part2.getCharge()!=0 ) { 
                                if(pass==1)
                                    if(this.checkVertex(part1,part2)==false) {
                                        continue;
                                    }
                                if(pass>1)
                                    if(this.checkSequentialVertex(part1,part2)==false) {
                                        continue;
                                    }
                                if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)==false) continue;
                            }
                            
                                
                            if(part1.getCharge()!=0 && part2.getCharge()==0) { 
                                double[] p1 = part1.SwimTo(0,part2.getVx(), part2.getVy(), part2.getVz(), 
                                        part2.getPx(), part2.getPy(), part2.getPz(), swim);
                                double vx1 = p1[0];
                                double vy1 = p1[1];
                                double vz1 = p1[2];
                                double px1 = p1[3];
                                double py1 = p1[4];
                                double pz1 = p1[5];
                                part1.setVx(vx1);
                                part1.setVy(vy1);
                                part1.setVz(vz1);
                                part1.setPx(px1);
                                part1.setPy(py1);
                                part1.setPz(pz1);
                                if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)==false) continue;

                            } else if(part2.getCharge()!=0 && part1.getCharge()==0) { 
                                double[] p2 = part2.SwimTo(0,part1.getVx(), part1.getVy(), part1.getVz(), 
                                        part1.getPx(), part1.getPy(), part1.getPz(), swim);
                                double vx2 = p2[0];
                                double vy2 = p2[1];
                                double vz2 = p2[2];
                                double px2 = p2[3];
                                double py2 = p2[4];
                                double pz2 = p2[5];
                                part2.setVx(vx2);
                                part2.setVy(vy2);
                                part2.setVz(vz2);
                                part2.setPx(px2);
                                part2.setPy(py2);
                                part2.setPz(pz2);
                               if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)==false) continue;
                            } else if(part2.getCharge()==0 && part1.getCharge()==0) { 
                                //System.out.println("PID1 "+part1.getPid()+" PID2 "+part2.getPid());
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
                                double vx1 = parExtrap.origin().x();
                                double vy1 = parExtrap.origin().y();
                                double vz1 = parExtrap.origin().z();
                                double vx2 = parExtrap.end().x();
                                double vy2 = parExtrap.end().y();
                                double vz2 = parExtrap.end().z();
                                part1.setVx(vx1);
                                part1.setVy(vy1);
                                part1.setVz(vz1);
                                part2.setVx(vx2);
                                part2.setVy(vy2);
                                part2.setVz(vz2);
                                if(part.combine(part1, part2, parPID, loMassCut, hiMassCut)==false) continue;

                            }
                            part1.isUsed = true;
                            part2.isUsed = true;
                            
                            //System.out.println("added part "+part.toString()+" with" );
                            //System.out.println("daugh1 "+part.getDaughters().get(0).toString());
                            //System.out.println("daugh2 "+part.getDaughters().get(1).toString());
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
        listsByPID.clear(); //System.out.println("IN SORTING for "+dau1PID+" "+dau2PID);
        for(Particle par : daughters) { //System.out.println(par.getPid());
            if(par.getPid()==dau1PID || par.getPid()==dau2PID || par.getPid()==dau3PID) {
                if(listsByPID.containsKey(par.getPid())) {
                    listsByPID.get(par.getPid()).add(par); //System.out.println("added to list "+par.getPid());
                } else {
                    listsByPID.put(par.getPid(), new ArrayList<Particle>());
                    listsByPID.get(par.getPid()).add(par);  //System.out.println("added to list "+par.getPid());
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
    private DoubleSwim _ds;
    private boolean checkSequentialVertex(Particle p1, Particle p2) {
        boolean pass = false;
        
         _ds.init(p1.getVx(), p1.getVy(), p1.getVz(),
            p1.getPx(), p1.getPy(), p1.getPz(), p1.getCharge(),
            p2.getVx(), p2.getVy(), p2.getVz(),
            p2.getPx(), p2.getPy(), p2.getPz(), p2.getCharge());

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
    
    private boolean checkVertex(Particle p1, Particle p2) {
        boolean pass = false;
        if(getVertBank()==null) return pass;
        
        if(getVertBank()!=null) {
            int nrows2 = getVertBank().rows();
            for(int loop2 = 0; loop2 < nrows2; loop2++){
                //double r = 999;
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
                    
                    //r =  (double) getVertBank().getFloat("r", loop2);
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
                    
                    pass=true;
                    return pass;
                }
                
            }
            
        }
       
        return pass;
    }

    private boolean overlaps(Particle part01, Particle part02) {
        boolean overlaps = false;
        
        List<Particle> part1 = new ArrayList<>();
        List<Particle> part2 = new ArrayList<>();
        part1.add(part01);
        part2.add(part02);
        System.out.println("ov"+part01.toString());
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

    
}
