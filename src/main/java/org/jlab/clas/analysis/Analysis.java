package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlab.clas.analysis.event.Reader;
import org.jlab.clas.analysis.event.Writer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.clas.swimtools.Swim;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.calib.utils.RCDBConstants;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
/**
 * Service to return reconstructed vertices from EB tracks
 *
 * @author ziegler
 *
 */


public class Analysis extends ReconstructionEngine {

    public Analysis(String name) {
        super(name, "ziegler", "1.0");
    }
    
    String FieldsConfig = "";
    private int Run = -1;
    
    //e.g. 3122:2212:-211:0:1.0:1.5;-3122:-2212:211:0:1.0:1.5
    private List<Integer> parent;
    private List<Integer> dau1;
    private List<Integer> dau2;
    private List<Integer> dau3;
    private List<Double> lowBound;
    private List<Double> highBound;
    private List<Integer> noBeamConst= new ArrayList<>();
    RCDBConstants rcdb = null;
    ConstantsManager conman = new ConstantsManager();
    
    public void setDecays(String decays) {
        parent      = new ArrayList<>();
        dau1        = new ArrayList<>();
        dau2        = new ArrayList<>();
        dau3        = new ArrayList<>();
        lowBound    = new ArrayList<>();
        highBound   = new ArrayList<>();
        
        if(decays!=null) {
            String[] chain = decays.split(";");
            for(int i =0; i< chain.length; i++) {
                String[] decprod = chain[i].split(":");
                System.out.println(decprod[0]+" --> "+decprod[1]+" + "+decprod[2]);
                parent.add(Integer.parseInt(decprod[0]));
                dau1.add(Integer.parseInt(decprod[1]));
                dau2.add(Integer.parseInt(decprod[2]));
                dau3.add(Integer.parseInt(decprod[3]));
                lowBound.add(Double.parseDouble(decprod[4]));
                highBound.add(Double.parseDouble(decprod[5]));
            }
        } else {
            System.out.println("!!!!!!!!!!!!!DECAYS ARE NULL!!!!!!!!!!!!!!!!!!!");
        }
    }

    public void updateWithNoBeamSpotPars(String nbspids) {
        noBeamConst = new ArrayList<>();
        
        if(nbspids!=null) {
            String[] chain = nbspids.split(":");
            for(int i =0; i< chain.length; i++) {
                System.out.println("use CVT no Beam Spot constraint for PID "+chain[i]);
                noBeamConst.add(Integer.parseInt(chain[i]));
            }
        } 
    }
    
    public int getRun() {
        return Run;
    }

    public void setRun(int run) {
        Run = run;
    }

    public String getFieldsConfig() {
        return FieldsConfig;
    }

    public void setFieldsConfig(String fieldsConfig) {
        FieldsConfig = fieldsConfig;
    }
    
    public static double beamE = 10.6;
    private int pass =1;
    @Override
    public boolean processDataEvent(DataEvent event) {
        Reader reader ;
        Writer writer ;
        
        this.FieldsConfig = this.getFieldsConfig();
        if (event.hasBank("RUN::config") == false) {
            System.err.println("RUN CONDITIONS NOT READ!");
            return false;
        }
        
        reader = new Reader();
        writer = new Writer();
        
        int newRun = event.getBank("RUN::config").getInt("run", 0); 
        if (Run != newRun)  {
            this.setRun(newRun);
            if (event.getBank("RUN::config").getInt("run",0) >= 100) {
                rcdb = conman.getRcdbConstants(event.getBank("RUN::config").getInt("run",0));
                beamE=rcdb.getDouble("beam_energy");
            }
        } 
        //System.out.println("EVENT "+event.getBank("RUN::config").getInt("event", 0));
        this.Run = this.getRun();
        
        
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        DataBank recParts = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("ANAL::Particle"))         recParts   = event.getBank("ANAL::Particle");
        if(recRun == null) return true;
        if(recParts == null) {
            return true;
        } 
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        Swim swim = new Swim();
        reader.readBanks(event, pass, noBeamConst);
       
        List<Particle> allparts = new ArrayList<>();
        
        reader.getDaus().forEach((key,value) -> allparts.add(value));
        //System.out.println("CHECK DAUS "+allparts.size());
        List<Particle> parts = new ArrayList<>();
        for(int k = 0; k < this.parent.size(); k++) {
            Decay dec = new Decay(this.parent.get(k), 
                    this.dau1.get(k), 
                    this.dau2.get(k), 
                    this.dau3.get(k), 
                    this.lowBound.get(k), 
                    this.highBound.get(k), 
                    allparts, swim,pass);
            if(dec!=null) {
                if(dec.getParticles()!=null) {
                    parts.addAll(dec.getParticles()); 
                }
            }
        } 
        DataBank bank =null;
        if(parts.isEmpty()) {
            if(event.hasBank("ANAL::Particle")) event.removeBank("ANAL::Particle");
        } else {
            bank = Writer.fillBank(event, parts, "ANAL::Particle", pass);
        }
        if(bank!=null) {
            event.appendBanks(bank);
            //bank.show();
            //System.out.println("================================================");
            //System.out.println(event.getBank("ANALREC::Particle").getInt("event", 0));
            //System.out.println("================================================");
            //event.appendBank(bank);
        }
        //event.show();
        return true;
   }

    @Override
    public boolean init() {
       
        this.registerBanks();
        this.loadConfiguration();
        this.printConfiguration();
        return true;
    }
    private void registerBanks() {
        
        super.registerOutputBank("ANALREC::Particle");
    } 
    
    public static void main(String[] args) {
    }

    
    public void loadConfiguration() {            
        
        // general (pass-independent) settings
        String dec = "decays";
        if (this.getEngineConfigString(dec)!=null) 
            this.setDecays(this.getEngineConfigString(dec));
        if (this.getEngineConfigString("pass")!=null) {
            pass=Integer.parseInt(this.getEngineConfigString("pass"));
        }
        String nbs ="nobeamspot";
        if (this.getEngineConfigString(nbs)!=null) 
            this.updateWithNoBeamSpotPars(this.getEngineConfigString(nbs));
        Constants.Load();
    }
    
    
    public void printConfiguration() {            
        for(int i =0; i<this.parent.size(); i++) {
            System.out.println("["+this.getName()+"] Reconstructing "+
                    this.parent.get(i)+"-->"+this.dau1.get(i)+":"+this.dau2.get(i));   
        }
        
    }

    
    
}
