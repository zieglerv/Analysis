package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import org.jlab.clas.analysis.event.Reader;
import org.jlab.clas.analysis.event.Writer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.clas.swimtools.Swim;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.calib.utils.RCDBConstants;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
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
    int passn =0;
    int counter =0;
    //e.g. 3122:2212:-211:0:1.0:1.5;-3122:-2212:211:0:1.0:1.5
    private List<Integer> parent;
    private List<Integer> dau1;
    private List<Integer> dau2;
    private List<Integer> dau3;
    private List<Double> lowBound;
    private List<Double> highBound;
    private double lommass =-999.0;
    private double himmass =999.0;
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
                System.out.println("LOOKING FOR DECAY "+decprod[0]+" --> "+decprod[1]+" + "+decprod[2]);
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
 
    public  double beamE = 10.6;
    public  Reader.Target target = null;
    private int pass =1;
    @Override
    public boolean processDataEvent(DataEvent event) {
        Reader reader ;
        this.FieldsConfig = this.getFieldsConfig();
        if (event.hasBank("RUN::config") == false) {
            event.show();
            //System.err.println("RUN CONDITIONS NOT READ!");
            return false;
        }
        counter++;
        reader = new Reader();
        
        int newRun = event.getBank("RUN::config").getInt("run", 0); 
        if (Run != newRun)  {
            this.setRun(newRun); 
            if (event.getBank("RUN::config").getInt("run",0) >= 100) {
                rcdb = conman.getRcdbConstants(event.getBank("RUN::config").getInt("run",0));
                beamE=rcdb.getDouble("beam_energy")/1000.0;
                String tar = rcdb.getString("target");
                target = Reader.Target.valueOf(tar);
            }
        } 
        //System.out.println("pass"+pass+" EVENT "+event.getBank("RUN::config").getInt("event", 0) +" beamE "+beamE+" tar "+target.name());
        this.Run = this.getRun();
        
        
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recBankFT = null;
        DataBank recEvenEB = null;
        DataBank recParts = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("RECFT::Particle"))        recBankFT   = event.getBank("RECFT::Particle");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("ANAL::Particle"))         recParts   = event.getBank("ANAL::Particle");
        if(recRun == null) return true;
        if(pass>1 && recParts == null) {
            return true;
        } 
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        Swim swim = new Swim();
        reader.readBanks(event, pass);
       
        List<Particle> allparts = new ArrayList<>();
        
        reader.getDaus().forEach((key,value) -> allparts.add(value));
        List<Particle> parts = new ArrayList<>();
        for(int k = 0; k < this.parent.size(); k++) {
           
            Decay dec = new Decay(this.parent.get(k), 
                    this.dau1.get(k), 
                    this.dau2.get(k), 
                    this.dau3.get(k), 
                    this.lowBound.get(k), 
                    this.highBound.get(k), 
                    this.lommass, this.himmass, beamE, target,
                    allparts, reader.getElectron(), swim,pass);
            if(dec!=null) {
                if(dec.getParticles()!=null) {
                    parts.addAll(dec.getParticles()); 
                }
            }
        } 
        if(pass>1 && parts.isEmpty()) {
            return true;
        }
        if(pass==1 && !parts.isEmpty()) {
            DataBank evBank = event.createBank("ANAL::Event", 1);
            evBank.setFloat("tarmass",0, (float) target.getMass());
            evBank.setFloat("beamenergy",0, (float) beamE);
            event.appendBank(evBank);
        }
//        if(pass>1) {
//            for(Integer i : reader.getDaus().keySet()) {
//                if(reader.getDaus().get(i).getIdx()>99) {
//                    for(int k = 0; k < this.parent.size(); k++) {
//                        if(reader.getDaus().get(i).getPid()!=(int)this.dau1.get(k) && 
//                                reader.getDaus().get(i).getPid()!=(int)this.dau2.get(k) && 
//                                reader.getDaus().get(i).getPid()!=(int)this.dau3.get(k) ) {
//                            reader.getDaus().get(i).keepSameIdx=true;
//                            parts.add(reader.getDaus().get(i));
//                        }
//                    }
//                }
//            }
//        }
        DataBank bank =null;
        if(parts.isEmpty()) {
            if(event.hasBank("ANAL::Particle")) event.removeBank("ANAL::Particle");
        } else {
            bank = Writer.fillBank(event, parts, reader.getElectron(), "ANAL::Particle", pass);
        }
        if(bank!=null) {
            event.appendBanks(bank);
            //bank.show();
            passn++;
            //System.out.println("================================================");
            //System.out.println("PASS"+pass+" total so far "+counter+" passed "+passn);
            //System.out.println("================================================");
            //event.appendBank(bank);
        } else {
            ((HipoDataEvent) event).getHipoEvent().reset();
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
        
        super.registerOutputBank("ANAL::Particle");
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
        if (this.getEngineConfigString("lommass")!=null) {
            lommass=Double.parseDouble(this.getEngineConfigString("lommass"));
        }
        if (this.getEngineConfigString("himmass")!=null) {
            himmass=Double.parseDouble(this.getEngineConfigString("himmass"));
        }
        Constants.Load();
    }
    
    
    public void printConfiguration() {            
        for(int i =0; i<this.parent.size(); i++) {
            System.out.println("["+this.getName()+"] Reconstructing "+
                    this.parent.get(i)+"-->"+this.dau1.get(i)+":"+this.dau2.get(i));   
        }
        
    }

    
    
}
