import java.util.Arrays;

public class AminoAcidLL{


  //Class variables

  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  //Variables I created

  static AminoAcidLL head;
  static AminoAcidLL tail;

  //Temporary node used when sorting
  static AminoAcidLL sorted;

  //Count of how many nodes there are
  static int nodeCount = 0;

  //Counter used for several methods;
  static int counter = 0;

  //Static vars used for methods
  static char[] aminoAcids;
  static int[] aminoAcidsCounts;

  //Constructors

  AminoAcidLL(){

    next = null;

  }

  AminoAcidLL(String inCodon){

    //Determines what amino acid this codon corresponds to.
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);

    //Generates a list of codons that correspond to this aminoacid.
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);

    //Counts for how many times each possible codon appears in the RNA sequence.
    counts =  new int[codons.length];
    Arrays.fill(counts, 0);

    //Next aminoacid is null
    next = null;


  }


  //Methods

  /*Creates a new node, with a given amino acid/codon,
    pairs and increments the codon counter for that codon.
    NOTE: Does not check for repeats!! */
  private static void fillCounts(String codon){

    //Find the aminoAcid
    char aminoAcid = AminoAcidResources.getAminoAcidFromCodon(codon);

    //Iterate through LL
    AminoAcidLL temp = head;
    while (temp != null){

      //If we find aminoAcid(node) corresponding to this codon
      if(temp.aminoAcid == aminoAcid){

        //Iterate through codon list and find match
        for(int i = 0; i < temp.counts.length; i++){

          //Once match found => update count
          if(temp.codons[i].equals(codon)){

            temp.counts[i]++;

          }
        }
      }

      temp = temp.next;

    }
  }

  //Helper method that inserts new node at the end of LL
  private static void insertNewAminoAcid( AminoAcidLL newNode){

    //If list empty
    if(head == null) {
      head = newNode;
      tail = newNode;
      nodeCount++;
    }

    //If list isn't empty
    else {
      tail.next = newNode;
      tail = newNode;
      nodeCount++;
    }

  }

  /* Recursive method that increments the count for a specific codon:
     If it should be at this node, increments it and stops,
     if not passes the task to the next node.
     If there is no next node, add a new node to the list that would contain the codon.*/
  private static void addCodon(String codon){

    //Find what aminoacid this codon corresponds to
    char aminoAcid = AminoAcidResources.getAminoAcidFromCodon(codon);

    //If list is empty => insert aminoacid in list and fill count
    if(head == null) {
      insertNewAminoAcid(new AminoAcidLL(codon));
      fillCounts(codon);
    }

    //If list isn't empty => check for repeats
    else {

      //Iterate through list
      AminoAcidLL temp = head;
      while (temp.next != null) {

        //If repeat fill codon count
        if (temp.aminoAcid == aminoAcid) {
          fillCounts(codon);
          return;
        }

        temp = temp.next;

      }

      //If no repeat then create a new node
      insertNewAminoAcid(new AminoAcidLL(codon));
      fillCounts(codon);

    }
  }

  // Static method for generating a linked list from an RNA sequence
  public static AminoAcidLL createFromRNASequence(String inSequence){

    //Iterate through sequence, get codons and pass them to the add codon method until done with sequence
    for(int i = 0; i < inSequence.length() - 2; i+=3){

      //Codon length 3
      String codon = inSequence.substring(i, i + 3);

      //Stop codons make the method stop
      if(codon.equals("UGA") || codon.equals("UAG") || codon.equals("UAA")){
        return head;
      }

      //Add codon to list
      addCodon(codon);

    }

    // Return the first node
    return head;

  }

  //Helper method to print nodes and their contents
  public static void printNodes(){

    AminoAcidLL temp = head;

    while(temp != null){

      System.out.println("\nAmino acid: " + temp.aminoAcid);
      System.out.println("Codons and count for each:");
      for(int i = 0; i< temp.codons.length;i++){

        System.out.println(temp.codons[i] + " is used " + temp.counts[i] + " times.");

      }

      temp = temp.next;

    }

  }

  /* Shortcut to find the total number igt of instances of this amino acid */
  private int totalCount(){

    int count = 0;

    for(int i = 0; i < this.counts.length; i++){

      count +=counts[i];

    }

    return count;
  }

  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }

  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(inList.counts[i] - counts[i]);
    }
    return diff;
  }

  public int aminoAcidCompare(AminoAcidLL inList){

    //Return the difference
    return this.totalDiff(inList);

  }

  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){

    return codonDiff(inList);

  }


  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){

    //Give char[] it's length on first run
    if(counter == 0) {
      aminoAcids = new char[nodeCount];
    }

    //If last node, base case
    if(this.next == null){
      aminoAcids[counter] = this.aminoAcid;
      counter = 0;
      return aminoAcids;

    }

    //Add amino acid to array
    aminoAcids[counter] = this.aminoAcid;
    counter++;

    //Recursive call
    return next.aminoAcidList();

  }

  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){

    //Give int[] it's length on first run
    if(counter == 0) {
      aminoAcidsCounts = new int[nodeCount];
    }

    //If last node, base case
    if(this.next == null){
      aminoAcidsCounts[counter] = this.totalCount();
      counter = 0;
      return aminoAcidsCounts;

    }

    //Add count to array
    aminoAcidsCounts[counter] = this.totalCount();
    counter++;

    //Recursive call
    return next.aminoAcidCounts();

  }

  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){

    //If we reach end then list must be sorted
    if(next == null)
      return true;

    //If we find case where amino acids are out of order then we return false
    if(aminoAcid > next.aminoAcid) {
      return false;
    }

    //Recursive call
    return next.isSorted();

  }

  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){

    //Create temp node
    AminoAcidLL current = inList;

    //Iterate through list
    while (current != null){

      //Save the next node
      AminoAcidLL nextNode = current.next;

      //Insert node where it belongs
      sortedInsert(current);

      //Next node
      current = nextNode;

    }

    head = sorted;
    return head;

  }

  //helper method for sort
  private static void sortedInsert(AminoAcidLL newnode){

    //If first node to be sorted, list will be sorted in new LL with head being called 'sorted'
    if(sorted == null || sorted.aminoAcid >= newnode.aminoAcid){

      newnode.next = sorted;
      sorted = newnode;

    }

    else{

      AminoAcidLL current = sorted;

      while(current.next != null && current.next.aminoAcid < newnode.aminoAcid){

        current = current.next;

      }

      newnode.next = current.next;
      current.next = newnode;

    }
  }

}