package taxon;

public class Taxon {

    public static int count = 0;

    public int id;
    public String label;
    
    
    public Taxon(int i, String lb){
        id = i;
        label = lb;
    }

    public Taxon(String lb){
        id = count++;
        label = lb;
    }

    public static String normalizeLabel(String label) {
        int separatorIndex = label.indexOf('_');
        if (separatorIndex < 0) {
            return label;
        }
        return label.substring(0, separatorIndex);
    }

    @Override
    public String toString(){
        return label;
    }

    
} 
