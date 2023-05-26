/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis.event;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlab.clas.analysis.Decay;
import org.jlab.clas.analysis.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
/**
 *
 * @author ziegler
 */
public class Reader {
    private Map<Integer, Particle> parts = new HashMap<>();
    private Map<Integer, Particle> daus = new HashMap<>();
    public static boolean useMCTruth = true;
    public boolean iseDetected = false;
    private final boolean updateWithUTrack = true;
   
    public void readDataBanksold(DataEvent event) {
        getParts().clear();
        getDaus().clear();
        List<Particle> parts = this.getPartsFromRec(event);
        if(parts!=null) {
            for(Particle part : parts) {
                getDaus().put(part.getIdx(), part);
            }
        }
    }
    
    public List<Particle> getPartsFromRec(DataEvent event) {
        iseDetected = false;
        DataBank runConf = null; 
        DataBank recBankEB = null;
        DataBank recDeteEB = null; 
        DataBank vertBankEB = null;
        DataBank trkBankEB = null;
        DataBank utrkBankEB = null;
        List<Particle> parts = new ArrayList<>();
        int ev = 0;
        if(event.hasBank("REC::Particle")) recBankEB = event.getBank("REC::Particle");
        if(event.hasBank("REC::Track")) trkBankEB = event.getBank("REC::Track");
        if(event.hasBank("REC::UTrack")) utrkBankEB = event.getBank("REC::UTrack");
        if(event.hasBank("REC::Calorimeter")) recDeteEB = event.getBank("REC::Calorimeter");
        if(event.hasBank("REC::VertDoca")) vertBankEB = event.getBank("REC::VertDoca");
        Decay.setVertBank(vertBankEB);        
        
        if(event.hasBank("RUN::config")) {   
            runConf   = event.getBank("RUN::config"); 
            ev = runConf.getInt("event", 0);
        }
        Map<Integer, double[]> uTrkMap = new HashMap<>();
        Map<Integer, double[]> pTrkMap = new HashMap<>();
        if(utrkBankEB!=null) {
            int nrows2 = utrkBankEB.rows();
            for(int loop = 0; loop < nrows2; loop++){
                int uindex = utrkBankEB.getInt("index", loop);
                double px = utrkBankEB.getFloat("px", loop);
                double py = utrkBankEB.getFloat("py", loop);
                double pz = utrkBankEB.getFloat("pz", loop);
                double vx = utrkBankEB.getFloat("vx", loop);
                double vy = utrkBankEB.getFloat("vy", loop);
                double vz = utrkBankEB.getFloat("vz", loop);
                double[] t = new double[]{px,py,pz,vx,vy,vz};
                uTrkMap.put(uindex, t);
            }
        }
        if(recBankEB!=null) {
            if(recBankEB.getInt("pid", 0)==11)
                    this.iseDetected = true;
            int nrows = recBankEB.rows();
            for(int loop = 0; loop < nrows; loop++){
                double px = recBankEB.getFloat("px", loop);
                double py = recBankEB.getFloat("py", loop);
                double pz = recBankEB.getFloat("pz", loop);
                double vx = recBankEB.getFloat("vx", loop);
                double vy = recBankEB.getFloat("vy", loop);
                double vz = recBankEB.getFloat("vz", loop);
                double[] t = new double[]{px,py,pz,vx,vy,vz};
                pTrkMap.put(loop, t);
            }
        }
        if(trkBankEB!=null) {
            int nrows2 = trkBankEB.rows();
            for(int loop = 0; loop < nrows2; loop++){
                int detector = trkBankEB.getInt("detector", loop);
                if(detector!=5) 
                    continue;
                int index = trkBankEB.getInt("index", loop);
                int pindex = trkBankEB.getInt("pindex", loop);
                if(uTrkMap.containsKey(index) && pTrkMap.containsKey(pindex)) {
                    pTrkMap.put(pindex, uTrkMap.get(index));
                }
            }
        }
        if(recBankEB!=null) {
            int nrows = recBankEB.rows();
            for(int loop = 0; loop < nrows; loop++){
                int pidCode = recBankEB.getInt("pid", loop);
                if(pidCode==0) continue;
                double px = recBankEB.getFloat("px", loop);
                double py = recBankEB.getFloat("py", loop);
                double pz = recBankEB.getFloat("pz", loop);
                double vx = recBankEB.getFloat("vx", loop);
                double vy = recBankEB.getFloat("vy", loop);
                double vz = recBankEB.getFloat("vz", loop);
                int det = -1;
                
                if((int) (Math.abs(recBankEB.getInt("status", loop))/1000)==4) 
                    det = 0;
                if((int) (Math.abs(recBankEB.getInt("status", loop))/1000)==2) 
                    det = 1;
                if(this.updateWithUTrack) {
                    double[] t = pTrkMap.get(loop);
                    px = t[0];
                    py = t[1];
                    pz = t[2];
                    vx = t[3];
                    vy = t[4];
                    vz = t[5];
                } 
                org.jlab.clas.physics.Particle recParticle = new org.jlab.clas.physics.Particle(
                                                pidCode,
                                                px,
                                                py,
                                                pz,
                                                vx,
                                                vy,
                                                vz);
                if(pidCode==22 && recDeteEB!=null) {
                    double energy1=0;
                    double energy4=0;
                    double energy7=0;
                    int    sector =0;
                    int  detector =0;
                    for(int j=0; j<recDeteEB.rows(); j++) {
                        if(recDeteEB.getShort("pindex",j)==loop && recDeteEB.getByte("detector",j)==DetectorType.ECAL.getDetectorId()) {
                            detector = recDeteEB.getByte("detector",j);
                            if(energy1 >= 0 && recDeteEB.getByte("layer",j) == 1) {
                                energy1 += recDeteEB.getFloat("energy",j);
                                sector = recDeteEB.getByte("sector",j);
                            }
                            if(energy4 >= 0 && recDeteEB.getByte("layer",j) == 4) energy4 += recDeteEB.getFloat("energy",j);
                            if(energy7 >= 0 && recDeteEB.getByte("layer",j) == 7) energy7 += recDeteEB.getFloat("energy",j);
                        }
                    }
                    recParticle.setProperty("energy1",energy1);
                    recParticle.setProperty("energy4",energy4);
                    recParticle.setProperty("energy7",energy7);
                    recParticle.setProperty("sector",sector*1.0);
                    recParticle.setProperty("detector",detector*1.0);
                    if(recParticle.charge()==0 
                            && recParticle.getProperty("detector")==DetectorType.ECAL.getDetectorId()) {
                        if(energy1>0.05 && energy4>0.0) {
                            double energy=(energy1+energy4+energy7)/0.245;
                            recParticle.setProperty("energy", energy);
                            Particle part = new Particle(recParticle);
                            part.setIdx(loop+1);
                            part.setPid(pidCode);
                            part.setUvx(vx);
                            part.setUvy(vy);
                            part.setUvz(vz);
                            part.setUpx(px);
                            part.setUpy(py);
                            part.setUpz(pz);
                            part.setDet(det);
                            parts.add(part);
                        }
                    }
                }
                
                if(Math.abs(pidCode)==211 || Math.abs(pidCode)==321  || Math.abs(pidCode)==2212 || pidCode==11) {
                    
                    double beta = (double)recBankEB.getFloat("beta", loop);
                    if(beta>1.0) 
                        beta=1.0;
                    double calcBeta = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()
                            +recParticle.mass()*recParticle.mass());
                    double mass2   = Math.pow(recParticle.p()/beta, 2)-recParticle.p()*recParticle.p();
                    if(mass2<0) continue;
                    int status = (int) Math.abs(recBankEB.getShort("status", loop));
                    double chi2pid = (double) Math.abs(recBankEB.getFloat("chi2pid", loop));
                    recParticle.setProperty("status", (double) status);
                    recParticle.setProperty("chi2pid", (double) chi2pid);
                    recParticle.setProperty("beta", (double) beta);
                    recParticle.setProperty("calcbeta", (double) calcBeta);
                    recParticle.setProperty("mass", recParticle.mass()); 
                    
                    if(chi2pid<Double.POSITIVE_INFINITY ) {    
                        Particle part = new Particle(recParticle);
                        
                        if(pidCode==11 && Math.abs(status)<4000)  { 
                            Particle parte = new Particle(recParticle); 
                            parte.setIdx(loop+1);
                            parte.setPid(pidCode);
                            parte.setUmass(Math.sqrt(mass2)); 
                            parte.setUvx(vx);
                            parte.setUvy(vy);
                            parte.setUvz(vz);
                            parte.setUpx(px);
                            parte.setUpy(py);
                            parte.setUpz(pz);
                            parte.setDet(det);
                            parts.add(parte);
                        }
                        if(Math.abs(part.getCharge())==1 && 
                               ( Math.abs(pidCode)==2212 || Math.abs(pidCode)==211 || Math.abs(pidCode)==321 )  ) {
                            Particle parth = new Particle(recParticle); 
                            parth.setIdx(loop+1);
                            parth.setPid(pidCode);
                            parth.setUmass(Math.sqrt(mass2)); 
                            parth.setUvx(vx);
                            parth.setUvy(vy);
                            parth.setUvz(vz);
                            parth.setUpx(px);
                            parth.setUpy(py);
                            parth.setUpz(pz);
                            parth.setDet(det);
                            parts.add(parth); //System.out.println("READ "+pidCode);
                        }
                    }
                }
            }
        } 
        return parts;
    }
    public void readBanks(DataEvent de) {
        readAnalBanks(de);
        List<Particle> parts = this.getPartsFromRec(de);
        if(parts!=null) {
            for(Particle part : parts) {
                getDaus().put(part.getIdx(), part);
            }
        }
        getDaus().forEach((key,value)->System.out.println("check "+key+" "+ value.getPid()));
    }
    public void readAnalBanks(DataEvent de) {
        getParts().clear();
        getDaus().clear();
        DataBank partBank = null;
        if(de.hasBank("ANAL::Particle")) partBank = de.getBank("ANAL::Particle");
        if(partBank == null) 
            return;
        Map<Integer, Particle> partMap = new HashMap<>();
        for(int i=0; i<partBank.rows(); i++) {
            int ev = partBank.getInt("event", i);
            int idx = partBank.getShort("idx",i); System.out.println("read idx "+idx);
            int pid = partBank.getInt("pid",i);
            double e = partBank.getFloat("e",i);
            double px = partBank.getFloat("px",i);
            double py = partBank.getFloat("py",i);
            double pz = partBank.getFloat("pz",i);
            double vx = partBank.getFloat("vx",i);
            double vy = partBank.getFloat("vy",i);
            double vz = partBank.getFloat("vz",i);
            double r = partBank.getFloat("r",i);
            int charge = partBank.getByte("charge",i);
            double mass = partBank.getFloat("mass",i);
            double umass = partBank.getFloat("umass",i);
            int ndau = partBank.getByte("ndau",i);
            int dau1idx = partBank.getShort("dau1idx",i);
            int dau2idx = partBank.getShort("dau2idx",i);
            int dau3idx = partBank.getShort("dau3idx",i);
            
            Particle part = new Particle( idx,  pid,  e,  px,  py,  pz, 
             0,  0,  0,  0,
             vx,  vy,  vz, 
             charge, mass,
             ndau,  dau1idx,  dau2idx,  dau3idx);
            int det = (int) partBank.getByte("det",i);
            part.setIdx(idx);
            part.setDet(det);
            part.setUmass(umass);
            part.setRecMass(mass);
            part.setR(r);
            part.setE(e);
            part.setMassConstrE(e);
            partMap.put(idx, part);
            if(ndau==0) {
                getDaus().put(idx, part);
            }
            if(ndau>1) {
                getParts().put(idx, part); 
            }
        }
        for (Integer i1 : partMap.keySet()) {
            for (Integer i2 : partMap.keySet()) {
                if(partMap.get(i1).hasDaughters() ){
                    for(int j = 0; j<partMap.get(i1).getDaughters().size(); j++) {
                        if(partMap.get(i1).getDaughters().get(j).getIdx()==partMap.get(i2).getIdx()) {
                            partMap.get(i1).getDaughters().set(j, partMap.get(i2));
                        }
                    }
                }
            }
        }
        for (Integer i1 : partMap.keySet()) {
            for (Integer i2 : partMap.keySet()) {
                if(partMap.get(i1).hasDaughters() ){
                    for(int j = 0; j<partMap.get(i1).getDaughters().size(); j++) {
                        if(partMap.get(i1).getDaughters().get(j).getIdx()==partMap.get(i2).getIdx()) {
                            partMap.get(i1).getDaughters().set(j, partMap.get(i2));
                        }
                    }
                }
            }
        }
        getDaus().clear();
        partMap.forEach((key,value) -> getDaus().put(key, value));
        getDaus().forEach((key,value)->System.out.println("...daus "+key+" "+ value.toString()));
 
//        while(it.hasNext()){
//           Map.Entry me = (Map.Entry)it.next();
//           Particle hpart = (Particle) me.getValue(); System.out.println("READ "+hpart.getPid());
//           int nd = hpart.getDaughters().size(); System.out.println("READ sz "+nd);
//           int d1x = hpart.getDaughters().get(0).getIdx();
//           int d2x = hpart.getDaughters().get(1).getIdx();
//           int d3x = -1;
//           if(nd == 3) 
//               d3x = hpart.getDaughters().get(2).getIdx();
//           hpart.getDaughters().clear();
//           hpart.getDaughters().add(getDaus().get(d1x));
//           hpart.getDaughters().add(getDaus().get(d2x));
//           if(nd == 3) 
//               hpart.getDaughters().add(getDaus().get(d3x));
//           
//        }
//        getDaus().clear();
//        getParts().forEach((key,value)->System.out.println("pars" +key+" "+ value.toString()));
//        getDaus().forEach((key,value)->System.out.println("...daus "+key+" "+ value.toString()));
//        getParts().forEach((key,value) -> getDaus().put(key, value));
//        getDaus().forEach((key,value)->System.out.println("...daus "+key+" "+ value.toString()));
    }

    
    /**
     * @return the parts
     */
    public Map<Integer, Particle> getParts() {
        return parts;
    }

    /**
     * @param parts the parts to set
     */
    public void setParts(Map<Integer, Particle> parts) {
        this.parts = parts;
    }

    /**
     * @return the daus
     */
    public Map<Integer, Particle> getDaus() {
        return daus;
    }

    /**
     * @param daus the daus to set
     */
    public void setDaus(Map<Integer, Particle> daus) {
        this.daus = daus;
    }
    boolean passEvent = false;
    private int getPidCode(DataBank mcBank, DataBank recBankEB, int loop) {
        int pid = recBankEB.getInt("pid", loop);
        if(mcBank!=null) {
            for(int loopm = 0; loopm < mcBank.rows(); loopm++){
                if(mcBank.getInt("pid", loopm) == 3122)
                    passEvent = true;
            }
        }
//        if(mcBank!=null) {
//            double px = recBankEB.getFloat("px", loop);
//            double py = recBankEB.getFloat("py", loop);
//            double pz = recBankEB.getFloat("pz", loop);
//            double p = Math.sqrt(px*px+py*py+pz*pz);
//            double phi = Math.atan2(py,px);
//            double theta = Math.acos(pz/p);
//            for(int loopm = 0; loopm < mcBank.rows(); loopm++){
//                if(getDaus().get(loop+1)!=null) {
//                    double pxmc = mcBank.getFloat("px", loopm);
//                    double pymc = mcBank.getFloat("py", loopm);
//                    double pzmc = mcBank.getFloat("pz", loopm);
//                    double pmc = Math.sqrt(pxmc*pxmc+pymc*pymc+pzmc*pzmc);
//                    double phimc = Math.atan2(pymc,pxmc);
//                    double thetamc = Math.acos(pzmc/pmc);
//                    if(Math.abs(pmc-p)/pmc<5*0.05 
//                            && Math.abs(phimc-phi)<5*0.005
//                            && Math.abs(thetamc-theta)<5*0.0010) {
//                        pid = (int) mcBank.getInt("pid", loopm);
//                    }
//                }  
//            }
//        }
        return pid;
    }

    private org.jlab.clas.physics.Particle truthMatch(org.jlab.clas.physics.Particle recParticle, List<org.jlab.clas.physics.Particle> mcParts) {
        org.jlab.clas.physics.Particle matched = null;
        double angle0 = 999.0;
        for(int i = 0; i < mcParts.size(); i++) {
            if(mcParts.get(i).pid()!=recParticle.pid()) 
                continue;
            Vector3D recP = new Vector3D(recParticle.px(), recParticle.py(), recParticle.pz());
            Vector3D genP = new Vector3D(mcParts.get(i).px(), mcParts.get(i).py(), mcParts.get(i).pz());
            double angle = recP.asUnit().angle(genP.asUnit());
            if(angle<angle0) {
                angle0 = angle;
                matched = mcParts.get(i);
            }
        }
        return matched;
    }

    public void updateDausList(List<Particle> parts) {
        for(Particle part : parts) {
            if(!getDaus().containsKey(part.getIdx())) {
                getDaus().put(part.getIdx(), part);
            }
        }
        
    }

}
