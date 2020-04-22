import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;

class AminoAcidResourcesTester{

  //Tests that were included
  @Test
  public void allCodons(){
    char[] rna = {'A','C','U','G'};
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        for(int k=0; k<4;k++){
          String s = new String(new char[]{rna[i],rna[j],rna[k]});
          char aaOut = AminoAcidResources.getAminoAcidFromCodon(s);
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }

          aaOut = AminoAcidResources.getAminoAcidFromCodon(s.toLowerCase());
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }
        }
      }
    }

  }

  @Test
  public void allAAs(){

    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<aa.length; i++){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }

      codonList = AminoAcidResources.getCodonListForAminoAcid(Character.toLowerCase(aa[i]));
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }
    }
  }

  //My test cases, 5 assertions per test.
  @Test
  public void Test1(){

    //Create from RNA method
    AminoAcidLL.createFromRNASequence("CCGUUGGCACUGUUG");
    System.out.println("### Initial linked list ###");
    AminoAcidLL.printNodes();

    //Sort method
    AminoAcidLL.sort(AminoAcidLL.head);
    System.out.println("\n### Sorted linked list ###");
    AminoAcidLL.printNodes();

    //Is sorted? method
    boolean sorted = AminoAcidLL.head.isSorted();
    assertTrue(sorted);

    //Amino Acid list method
    char[] aminoAcidsList = AminoAcidLL.head.aminoAcidList();
    char[] expected = {'A', 'L', 'P'};
    assertArrayEquals(expected, aminoAcidsList);

    //Amino acid counts method
    int[] countsList = AminoAcidLL.head.aminoAcidCounts();
    int[] expected2 = {1 ,3 ,1};
    assertArrayEquals(expected2, countsList);

    //Amino acid compare method
    int comparedDifference = AminoAcidLL.head.aminoAcidCompare(AminoAcidLL.tail);
    int expectedDifference = 0;
    assertEquals(expectedDifference, comparedDifference);

    //Codon compare method
    int comparedDifference2 = AminoAcidLL.head.codonCompare(AminoAcidLL.head.next);
    int expectedDifference2 = 2;
    assertEquals(expectedDifference2, comparedDifference2);

  }

  @Test
  public void Test2(){

    //Create from RNA method
    AminoAcidLL.createFromRNASequence("UGGGCGUGCAAG");
    System.out.println("### Initial linked list ###");
    AminoAcidLL.printNodes();

    //Sort method
    AminoAcidLL.sort(AminoAcidLL.head);
    System.out.println("\n### Sorted linked list ###");
    AminoAcidLL.printNodes();

    //Is sorted? method
    boolean sorted = AminoAcidLL.head.isSorted();
    assertTrue(sorted);

    //Amino Acid list method
    char[] aminoAcidsList = AminoAcidLL.head.aminoAcidList();
    char[] expected = {'A', 'C', 'K', 'W'};
    assertArrayEquals(expected, aminoAcidsList);

    //Amino acid counts method
    int[] countsList = AminoAcidLL.head.aminoAcidCounts();
    int[] expected2 = {1 ,1 ,1, 1};
    assertArrayEquals(expected2, countsList);

    //Amino acid compare method
    int comparedDifference = AminoAcidLL.head.aminoAcidCompare(AminoAcidLL.tail);
    int expectedDifference = 0;
    assertEquals(expectedDifference, comparedDifference);

    //Codon compare method
    int comparedDifference2 = AminoAcidLL.head.codonCompare(AminoAcidLL.head);
    int expectedDifference2 = 0;
    assertEquals(expectedDifference2, comparedDifference2);

  }


}