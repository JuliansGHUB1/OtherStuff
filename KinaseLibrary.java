package org.panda.causalpath.resource;

import org.panda.resource.FileServer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class KinaseLibrary extends FileServer {

    /*
    A map which can be used to determine the probability that a kinase has some
    amino acid at some position
     */
    private Map<String, Map<Integer, Map<String, Double>>> probabilityMap;

    @Override
    public String[] getLocalFilenames() {
        return new String[]{"KinaseLibraryMatrices.txt"};
    }


    // Returns all kinases within the map
    public Set<String> kinaseSet(){
        return probabilityMap.keySet();
    }

    /* Public method for user to access the probability that if kinase phosphorlyated
     a given peptide, what are the odds the peptide has given aminoAcid at specified location
     */
    public double aminoAcidProbability(String kinase, int location, String aminoAcid){
        return probabilityMap.get(kinase).get(location).get(aminoAcid);
    }

    @Override
    public boolean load() throws IOException {
        // Get the path to the file
        Path p = (Paths.get(this.locateInBase(this.getLocalFilenames()[0])));

        // Read all the lines to begin processing
        List<String> lines = Files.readAllLines(p);

        // Process the first line in order to get column headers (do not hardcode these)
        String[] firstLineSplit = lines.get(0).split("\t");

        // This next part of the code will find the locations/amino acids from the column headers
        // We will store these in two arrays

        int[] locations = new int[lines.size()];
        String[] aminoAcids = new String[lines.size()];

        splitColumnHeaders(firstLineSplit, locations, aminoAcids);

        // Process row-by-row to fill out maps

         /* A 3 layered map, which can be used to access raw data regarding a kinase's
       preference for a specific amino acid at some location
        */
        Map<String, Map<Integer, Map<String, Double>>> kinaseLocationMap = new HashMap<>();

        probabilityMap = new HashMap<>();

        for(int row = 1; row < lines.size(); row++){
            // Split the current row or line by the delimiter tab
            String[] splitLine = lines.get(row).split("\t");

            // Determine the kinase for curr row, and put it into the kinase-Location map
            String currKinase = splitLine[0];

            kinaseLocationMap.put(currKinase, new HashMap<>());
            probabilityMap.put(currKinase, new HashMap<>());

            for(int cell = 1; cell < splitLine.length; cell++){

                // Use the column headers to determine the location/amino acid for current cell
                int currLocation = locations[cell];
                String currAminoAcid = aminoAcids[cell];
                double measurement = Double.parseDouble(splitLine[cell]);

                if(!kinaseLocationMap.get(currKinase).containsKey(currLocation)){
                    kinaseLocationMap.get(currKinase).put(currLocation, new HashMap<>());
                    probabilityMap.get(currKinase).put(currLocation, new HashMap<>());
                }

                kinaseLocationMap.get(currKinase).get(currLocation).put(currAminoAcid, measurement);
                probabilityMap.get(currKinase).get(currLocation).put(currAminoAcid, measurement);

            }

            /* Once we have processed a single row, we may go ahead and calculate the probability that
               an amino acid has occured at a given location
             */

            // Iterate through all locations for this kinase
            for(int location: kinaseLocationMap.get(currKinase).keySet()){
                // For each location for the given kinase, we will iterate through aa-raw measurement map
                Map<String, Double> aminoAcidValue = kinaseLocationMap.get(currKinase).get(location);

                // A variable to store the sum of the measurements for a given location
                double sum = 0;

                // Iterate through all amino-acid measurements for a given location and sum them
                for(double d: aminoAcidValue.values()){
                    sum = sum + d;
                }

                // For each amino acid, compute the probability by dividing corresponding value
                // by the sum of the measurements for that site
                for(String aa: aminoAcidValue.keySet()){

                    // Test
                    if(currKinase.equals("AURA") && location == -5){
                        System.out.println("sum for AURA is: " + sum);
                        if(sum - 183554785.61 < 0.05 && sum - 183554785.61 > -0.05){
                            System.out.println("Resonable margin of error");
                        }
                        else{
                            System.out.println(sum - 183554785.61);
                        }
                    }
                    else if(currKinase.equals("ATM") && location == 3){
                        System.out.println("sum for ATM is: " + sum);
                        if(sum - 137395297.7 < 0.05 && sum - 137395297.7 > -0.05){
                            System.out.println("Resonable margin of error");
                        }
                        else{
                            System.out.println(sum - 183554785.61);
                        }
                    }


                    // Delete later


                    double aaMeasurement = probabilityMap.get(currKinase).get(location).get(aa);
                    double aaProbability = aaMeasurement/sum;
                    probabilityMap.get(currKinase).get(location).put(aa, aaProbability);
                }


            }


        }

        return true;
    }


    public static void main(String[] args){

        // Load

        KinaseLibrary kN = new KinaseLibrary();
        try{
            kN.load();
        }
        catch(Exception e){

        }

        // Currently have some bad testing in here. Will flesh out and do it in a better way, just wanted to get some ideas.

        // Test #1: Testing whether splitColumnHeaders class works

        // Start with all the column headers
        String[] a = {"-5P", "-5G", "-5A", "-5C", "-5S", "-5T", "-5V", "-5I", "-5L", "-5M", "-5F", "-5Y", "-5W", "-5H", "-5K", "-5R", "-5Q", "-5N", "-5D", "-5E", "-5s", "-5t", "-5y", "-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N",
                "-4D", "-4E", "-4s", "-4t", "-4y", "-3P", "-3G", "-3A", "-3C", "-3S", "-3T", "-3V", "-3I", "-3L", "-3M", "-3F", "-3Y", "-3W", "-3H", "-3K", "-3R", "-3Q", "-3N", "-3D", "-3E", "-3s", "-3t", "-3y", "-2P", "-2G", "-2A", "-2C", "-2S", "-2T", "-2V", "-2I", "-2L", "-2M", "-2F", "-2Y", "-2W", "-2H", "-2K", "-2R", "-2Q", "-2N", "-2D", "-2E", "-2s", "-2t", "-2y", "-1P", "-1G", "-1A", "-1C",
                "-1S", "-1T", "-1V", "-1I", "-1L", "-1M", "-1F", "-1Y", "-1W", "-1H", "-1K", "-1R", "-1Q", "-1N", "-1D", "-1E", "-1s", "-1t", "-1y", "1P", "1G", "1A", "1C", "1S", "1T", "1V", "1I", "1L", "1M", "1F", "1Y", "1W", "1H", "1K", "1R", "1Q", "1N", "1D", "1E", "1s", "1t", "1y", "2P", "2G", "2A", "2C", "2S", "2T", "2V", "2I", "2L", "2M", "2F", "2Y", "2W", "2H", "2K", "2R", "2Q", "2N", "2D", "2E", "2s", "2t", "2y",
                "3P", "3G", "3A", "3C", "3S", "3T", "3V", "3I", "3L", "3M", "3F", "3Y", "3W", "3H", "3K", "3R", "3Q", "3N", "3D", "3E", "3s", "3t", "3y", "4P", "4G", "4A", "4C", "4S", "4T", "4V", "4I", "4L", "4M", "4F", "4Y", "4W", "4H", "4K", "4R", "4Q", "4N", "4D", "4E", "4s", "4t", "4y"};
        // Create 2 arrays for the locations/amino acids
        int[] locations = new int[a.length];
        String[] aminoAcids = new String[a.length];

        // Call the method
        kN.splitColumnHeaders(a, locations, aminoAcids);

        // Combine back locations and amino acids into s, and see if s is the same as a
        String[] s = new String[a.length];
        for(int i = 0; i < a.length; i++){
            s[i] = locations[i] + aminoAcids[i];
        }

        // Iterate through
        boolean test = true;
        for(int i = 0; i < a.length; i++){
            if(!a[i].equals(s[i])){
                test = false;
            }
        }

        System.out.println("splitColumnHeaders test passed: " + test);



        // Test #2 Testing whether the load correctly produces the matrices
        String[] columnHeader = {"-5P", "-5G", "-5A", "-5C", "-5S", "-5T", "-5V", "-5I", "-5L", "-5M", "-5F", "-5Y", "-5W", "-5H", "-5K", "-5R", "-5Q", "-5N", "-5D", "-5E", "-5s", "-5t", "-5y", "-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y", "-3P", "-3G", "-3A", "-3C", "-3S", "-3T", "-3V", "-3I", "-3L", "-3M", "-3F", "-3Y", "-3W", "-3H", "-3K", "-3R", "-3Q", "-3N", "-3D", "-3E", "-3s", "-3t", "-3y", "-2P", "-2G", "-2A", "-2C", "-2S", "-2T", "-2V", "-2I", "-2L", "-2M", "-2F", "-2Y", "-2W", "-2H", "-2K", "-2R", "-2Q", "-2N", "-2D", "-2E", "-2s", "-2t", "-2y", "-1P", "-1G", "-1A", "-1C", "-1S", "-1T", "-1V", "-1I", "-1L", "-1M", "-1F", "-1Y", "-1W", "-1H", "-1K", "-1R", "-1Q", "-1N", "-1D", "-1E", "-1s", "-1t", "-1y", "1P", "1G", "1A", "1C", "1S", "1T", "1V", "1I", "1L", "1M", "1F", "1Y", "1W", "1H", "1K", "1R", "1Q", "1N", "1D", "1E", "1s", "1t", "1y", "2P", "2G", "2A", "2C", "2S", "2T", "2V", "2I", "2L", "2M", "2F", "2Y", "2W", "2H", "2K", "2R", "2Q", "2N", "2D", "2E", "2s", "2t", "2y", "3P", "3G", "3A", "3C", "3S", "3T", "3V", "3I", "3L", "3M", "3F", "3Y", "3W", "3H", "3K", "3R", "3Q", "3N", "3D", "3E", "3s", "3t", "3y", "4P", "4G", "4A", "4C", "4S", "4T", "4V", "4I", "4L", "4M", "4F", "4Y", "4W", "4H", "4K", "4R", "4Q", "4N", "4D", "4E", "4s", "4t", "4y"
        };
        int[] locationStorage = new int[columnHeader.length];
        String[] aminoAcidStorage = new String[aminoAcids.length];

        kN.splitColumnHeaders(columnHeader, locationStorage, aminoAcidStorage);
        String currKinase = "AURA";
        String[] AURAMeasurements = {"5683377.86", "9080772.95", "8340840.22", "10883246.95", "8186154.23", "7555351.79", "5841037.6", "6726258.4", "5794366.62", "5560729.33", "6387081.93", "7418541.16", "6976863.18", "8319964.94", "10588197.24", "12875913.33", "6178303.19", "9242426.42", "7731860.37", "8130753.93", "8006580.23", "8006580.23", "10039583.51", "7958867.31", "10857162.09", "9473304.23", "11466629.34", "9472229.52", "7808415.78", "6145822.45", "6211134.39", "6735254.52", "7197377.31", "6378000.66", "6882198.26", "8495568.55", "8370676.86", "11899998.53", "15019256.73", "7990390.94", "9149076.87", "12280572.35", "8466437.59", "8550929.04", "8550929.04", "8461969.98", "7893630.45", "12341329.54", "5944089.31", "9439142.51", "9163620.4", "6872721.84", "3803083.24", "4114322.07", "3281708.07", "4837111.52", "5128332.52", "5133141.96", "6967519.69", "8615566.13", "14118447.67", "43189607.09", "5960485.23", "7663342.47", "5417073.9", "7601312.84", "7860535.17", "7860535.17", "8954818.68", "3594447.85", "3978731.61", "2608708.96", "3950832.04", "9450264.91", "4586496.58", "2023682.92", "2016602.06", "2266725.34", "2188762.8", "2530417.73", "2308419.4", "2883789.22", "4248220.34", "16373437.14", "146172558.6", "7567379.18", "4296821.91", "3407441.02", "3195616.4", "4660708.62", "4660708.62", "4432927.54", "6104437.77", "5618741.13", "5568670", "8329791.91", "8376059.15", "6398860.01", "4904681.74", "5198540.82", "10755585.68", "11285236.69", "7207505.51", "10853139.28", "4791017.71", "8897219.55", "5768673.36", "9428310.25", "8854185.69", "9383652.54", "11303201.82", "12748320.66", "4798936.23", "4798936.23", "18243735.87", "1570050.67", "4817906.3", "3530110.16", "13593499.82", "10018001.58", "6904104.97", "8300042.84", "12391701.34", "15242214.98", "13093979.88", "13596272.37", "8452920.14", "9417994.61", "4952508.91", "1794713.07", "3165114.16", "5469980.93", "6346712.55", "8523094.45", "5924121.94", "11965552.32", "11965552.32", "6221615.33", "6187683", "13923765.79", "7964465.97", "13763518.63", "9211406.25", "5923056.04", "6073157.6", "8137327.3", "6222667.55", "6597541.73", "7595807.16", "21996928.96", "6986969.28", "10348689.08", "2570453.2", "4039190.54", "4992070.17", "5148360.13", "7879384.3", "8182293.36", "10007316.82", "10007316.82", "12171087.94", "12164972.2", "13800089.39", "8459006.06", "14577163.74", "22444340.11", "7800360.93", "7234269.2", "8073261.69", "7576863.68", "7023019.46", "11680945.08", "12666320.57", "11212145.91", "12320805.5", "7874559.56", "7666434.04", "7864492.23", "7452652.89", "7673844.45", "7943308.09", "9963238.5", "9963238.5", "16243733.68", "8179334", "8900275.8", "7283974.75", "11543035.5", "9089075.19", "7825897.21", "5761667.26", "7859246.98", "5475055.46", "8466023.95", "6375030.33", "8001396.18", "10273142.04", "9946645.41", "10051110.92", "10434776.86", "8957963.67", "10614474.87", "7650278.9", "9278532.9", "12892094.54", "12892094.54", "12066131.31"};
/*
        boolean valid = true;
        for(int i = 0; i < locationStorage.length; i++){
            Double measurement = kN.kinaseLocationMap.get(currKinase).get(locationStorage[i]).get(aminoAcidStorage[i]);
            if(measurement != Double.parseDouble(AURAMeasurements[i])){
                valid = false;
            }
        }



        System.out.println("Result of Test2: " + valid);

        System.out.println(kN.probabilityMap.get("AURA").get(-5).get("P"));

        System.out.println(kN.probabilityMap.get("ATM").get(3).get("V"));

 */









    }



    private void splitColumnHeaders(String[] firstLine, int[] sites, String[] aminoAcids){
        for(int i = 0; i < firstLine.length; i++){
            int startIndex = 0;
            boolean foundDivider = false;
            while(startIndex < firstLine[i].length() && !foundDivider)
                if(firstLine[i].charAt(startIndex) == '-' || Character.isDigit(firstLine[i].charAt(startIndex))){
                    startIndex++;
                }
                else{
                    int siteValue = Integer.parseInt(firstLine[i].substring(0, startIndex));
                    String aminoAcidValue = firstLine[i].substring(startIndex);
                    sites[i] = siteValue;
                    aminoAcids[i] = aminoAcidValue;
                    foundDivider = true;
                }

        }
    }
}
