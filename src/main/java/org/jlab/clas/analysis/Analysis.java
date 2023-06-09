package org.jlab.clas.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlab.clas.analysis.event.Reader;
import org.jlab.clas.analysis.event.Writer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.service.vtx.VTXEngine;

/**
 * Service to return reconstructed vertices from EB tracks
 *
 * @author ziegler
 *
 */


public class Analysis extends ReconstructionEngine {

    public Analysis() {
        super("MYANAL", "ziegler", "1.0");
    }
    
    String FieldsConfig = "";
    private int Run = -1;
  
 

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
    
    public static boolean useVtxBank = false;
    
    
    boolean LOADED = false;
    
    private double zTarget;
    private double zLength;
    @Override
    public boolean processDataEvent(DataEvent event) {
        Reader reader ;
        Writer writer ;
        if(!event.hasBank("REC::VertDoca")) {
            VTXEngine vtxe = new VTXEngine();
            vtxe.init();
            vtxe.processDataEvent(event);
        }
        this.FieldsConfig = this.getFieldsConfig();
        if (event.hasBank("RUN::config") == false) {
            System.err.println("RUN CONDITIONS NOT READ!");
            return false;
        }
        reader = new Reader();
        writer = new Writer();
        
        int newRun = event.getBank("RUN::config").getInt("run", 0); 
        if (Run != newRun)  this.setRun(newRun); 
        //System.out.println("EVENT "+event.getBank("RUN::config").getInt("event", 0));
        this.Run = this.getRun();
        
        // Load target
        if(!LOADED) {
            ConstantProvider providerTG = GeometryFactory.getConstants(DetectorType.TARGET, this.getRun(), this.getConstantsManager().getVariation());
            this.zTarget = providerTG.getDouble("/geometry/target/position",0);
            this.zLength = providerTG.getDouble("/geometry/target/length",0);
            //System.out.println("TARGET   "+this.zTarget);
            LOADED = true;
        }
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(recRun == null) return true;
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        
        reader.readDataBanks(event, 0);
        //System.out.println("e detected "+reader.iseDetected);
        //DataBank runBank = ep.FillHeader(recRun);
        if(!reader.iseDetected)
            return false;
        List<Particle> allparts = new ArrayList<>();
        
        reader.getDaus().forEach((key,value) -> allparts.add(value));
        List<Particle> parts = new ArrayList<>();
        for(int k = 0; k < Constants.getInstance().parent.size(); k++) {
            Decay dec = new Decay(Constants.getInstance().parent.get(k), 
                    Constants.getInstance().dau1.get(k), 
                    Constants.getInstance().dau2.get(k), 
                    Constants.getInstance().dau3.get(k), 
                    Constants.getInstance().lowBound.get(k), 
                    Constants.getInstance().highBound.get(k), 
                    allparts);
            if(dec!=null) {
                if(dec.getParticles()!=null) {
                    parts.addAll(dec.getParticles());
                }
            }
        } 
        
        List<DataBank> banks = new ArrayList<>();
        banks.add(writer.fillBank(event, parts, ""));

        event.appendBanks(banks.toArray(new DataBank[0]));
        
        return true;
   }

    @Override
    public boolean init() {
        
        this.initConstantsTables();
        this.registerBanks();
        this.loadConfiguration();
        this.printConfiguration();
        return true;
    }
    private void registerBanks() {
        super.registerOutputBank("ANAL::Event");
        super.registerOutputBank("ANAL::Particle");
        super.registerOutputBank("ANAL::Daughter");
    } 
    
    public static void main(String[] args) {
    }

    
    public void loadConfiguration() {            
        
        // general (pass-independent) settings
        
        if (this.getEngineConfigString("decays")!=null) 
            Constants.getInstance().setDecays(this.getEngineConfigString("decays"));
        
        
    }


    public void initConstantsTables() {
        String[] tables = new String[]{
            "/geometry/beam/position",
            "/geometry/target"
        };
        requireConstants(Arrays.asList(tables));
        this.getConstantsManager().setVariation("default");
    }
    
    
    
    
    public void printConfiguration() {            
        for(int i =0; i<Constants.getInstance().parent.size(); i++) {
            System.out.println("["+this.getName()+"] Reconstructing "+
                    Constants.getInstance().parent.get(i)+"-->"+Constants.getInstance().dau1.get(i)+":"+Constants.getInstance().dau2.get(i));   
        }
        
    }

    
    
}
