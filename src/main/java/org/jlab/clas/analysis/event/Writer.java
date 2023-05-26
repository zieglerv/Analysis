/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis.event;

import java.util.List;
import org.jlab.clas.analysis.Particle;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author ziegler
 */
public class Writer {
    
    public static DataBank fillBank(DataEvent event, List<Particle> partlist, 
            String bankName, int pass) {
        if (partlist == null || partlist.isEmpty()) return null;
        
        int bankSize = 0;
        for(Particle p : partlist) {
            bankSize++;
            int nDau = p.getDaughters().size();
            bankSize+=nDau;
            for(Particle p2: p.getDaughters()) {
                bankSize+=p2.getDaughters().size();
            }
        }
        if(event.hasBank(bankName)) event.removeBank(bankName);
        DataBank partBank = event.createBank(bankName, bankSize);
        int ev = event.getBank("RUN::config").getInt("event", 0); 
        int i = -1;
        for (int ii = 0; ii < partlist.size(); ii++) {
            i++;
            partBank.setInt("event", i, ev);
            partBank.setShort("idx",i, (short) (pass*100+ii));
            Particle parent = partlist.get(ii);
            fillBankValues(partBank, i, parent);
            if(!parent.getDaughters().isEmpty()) {
                String ds;
                for(int j =0; j<partlist.get(ii).getDaughters().size(); j++) {
                    ds = "dau";
                    ds+=j+1;
                    ds+="idx";
                    partBank.setShort(ds,i, (short) parent.getDaughters().get(j).getIdx());
                    //System.out.println("Checkoutput "+partlist.get(ii).getDaughters().get(j).toString());
                }
                for(int j =0; j<parent.getDaughters().size(); j++) {
                    i++;
                    partBank.setInt("event", i, ev);
                    Particle dau1stgen = parent.getDaughters().get(j);
                    partBank.setShort("idx",i, (short) dau1stgen.getIdx());
                    fillBankValues(partBank, i, dau1stgen);
                    
                    for(int jj =0; jj<dau1stgen.getDaughters().size(); jj++) {
                        ds = "dau";
                        ds+=jj+1;
                        ds+="idx";
                        partBank.setShort(ds,i, (short) dau1stgen.getDaughters().get(jj).getIdx());
                        //System.out.println("Checkoutput "+partlist.get(ii).getDaughters().get(j).toString());
                    }
                    for(int jj =0; jj<dau1stgen.getDaughters().size(); jj++) {
                        i++;
                        partBank.setInt("event", i, ev);
                        Particle dau2ndgen = dau1stgen.getDaughters().get(jj);
                        partBank.setShort("idx",i, (short) dau2ndgen.getIdx());
                        fillBankValues(partBank, i, dau2ndgen);
                        
                    }
                }  
            }
        }
        
        //partBank.show();
        return partBank;

    }
    
    public static DataBank fillBankDebug(DataEvent event, List<Particle> partlist, String bankName) {
        if (partlist == null || partlist.isEmpty()) return null;
        System.out.println("FILLING BANK");
        int bankSize = 0;
        for(Particle p : partlist) {
            bankSize++;
            int nDau = p.getDaughters().size();
            bankSize+=nDau;
        }
        bankName = "ANAL::Particle";
        DataBank partBank = event.createBank(bankName, bankSize);
        int ev = event.getBank("RUN::config").getInt("event", 0); 
        int i = 0;
        for (int ii = 0; ii < partlist.size(); ii++) {
            
            partBank.setInt("event", i, ev);
            partBank.setShort("idx",i, (short) (100+ii));
            partBank.setInt("pid", i, partlist.get(ii).getPid());
            partBank.setFloat("e",i, (float) partlist.get(ii).getMassConstrE());
            partBank.setFloat("px",i, (float) partlist.get(ii).getPx());
            partBank.setFloat("py",i, (float) partlist.get(ii).getPy());
            partBank.setFloat("pz",i, (float) partlist.get(ii).getPz());
            partBank.setFloat("ecm",i, (float) partlist.get(ii).getECM());
            partBank.setFloat("pxcm",i, (float) partlist.get(ii).getPxcm());
            partBank.setFloat("pycm",i, (float) partlist.get(ii).getPycm());
            partBank.setFloat("pzcm",i, (float) partlist.get(ii).getPzcm());
            partBank.setFloat("vx",i, (float) partlist.get(ii).getVx());
            partBank.setFloat("vy",i, (float) partlist.get(ii).getVy());
            partBank.setFloat("vz",i, (float) partlist.get(ii).getVz());
            partBank.setFloat("r",i, (float) partlist.get(ii).getR());
            partBank.setByte("charge", i, (byte) partlist.get(ii).getCharge());
            partBank.setFloat("mass",i, (float) partlist.get(ii).getRecMass());
            partBank.setByte("ndau",i, (byte) partlist.get(ii).getDaughters().size());
            if(partlist.get(ii).getDaughters().size()>0) {
                String ds;
                for(int j =0; j<partlist.get(ii).getDaughters().size(); j++) {
                    ds = "dau";
                    ds+=j+1;
                    ds+="idx";
                    partBank.setShort(ds,i, (short) partlist.get(ii).getDaughters().get(j).getIdx());
                }
                for(int j =0; j<partlist.get(ii).getDaughters().size(); j++) {
                    i++;
                    partBank.setInt("event", i, ev);
                    partBank.setShort("idx",i, (short) partlist.get(ii).getDaughters().get(j).getIdx());
                    partBank.setInt("pid", i, partlist.get(ii).getDaughters().get(j).getPid());
                    partBank.setFloat("e",i, (float) partlist.get(ii).getDaughters().get(j).getMassConstrE());
                    partBank.setFloat("px",i, (float) partlist.get(ii).getDaughters().get(j).getPx());
                    partBank.setFloat("py",i, (float) partlist.get(ii).getDaughters().get(j).getPy());
                    partBank.setFloat("pz",i, (float) partlist.get(ii).getDaughters().get(j).getPz());
                    partBank.setFloat("ecm",i, (float) partlist.get(ii).getDaughters().get(j).getECM());
                    partBank.setFloat("pxcm",i, (float) partlist.get(ii).getDaughters().get(j).getPxcm());
                    partBank.setFloat("pycm",i, (float) partlist.get(ii).getDaughters().get(j).getPycm());
                    partBank.setFloat("pzcm",i, (float) partlist.get(ii).getDaughters().get(j).getPzcm());
                    partBank.setFloat("vx",i, (float) partlist.get(ii).getDaughters().get(j).getVx());
                    partBank.setFloat("vy",i, (float) partlist.get(ii).getDaughters().get(j).getVy());
                    partBank.setFloat("vz",i, (float) partlist.get(ii).getDaughters().get(j).getVz());
                    partBank.setByte("charge", i, (byte) partlist.get(ii).getDaughters().get(j).getCharge());
                    partBank.setFloat("mass",i, (float) partlist.get(ii).getDaughters().get(j).getRecMass());
                    partBank.setByte("ndau",i, (byte) partlist.get(ii).getDaughters().get(j).getDaughters().size());
                }
            }  
        }
        partBank.show();
        return partBank;

    }

    private static void fillBankValues(DataBank partBank, int i, Particle p) {
        partBank.setInt("pid", i, p.getPid());
        partBank.setFloat("e",i, (float) p.getMassConstrE());
        partBank.setFloat("px",i, (float) p.getPx());
        partBank.setFloat("py",i, (float) p.getPy());
        partBank.setFloat("pz",i, (float) p.getPz());
        partBank.setFloat("vx",i, (float) p.getVx());
        partBank.setFloat("vy",i, (float) p.getVy());
        partBank.setFloat("vz",i, (float) p.getVz());
        partBank.setFloat("r",i, (float) p.getR());
        partBank.setByte("charge", i, (byte) p.getCharge());
        partBank.setFloat("mass",i, (float) p.getRecMass());
        partBank.setFloat("umass",i, (float) p.getUmass());
        partBank.setByte("ndau",i, (byte) p.getDaughters().size());
        partBank.setByte("det", i, (byte) p.getDet());
    }
}
