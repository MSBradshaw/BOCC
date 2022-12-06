// this file is < 500 kb, really small not a problem to force someone's browswer to download and process ;)
var url = 'https://raw.githubusercontent.com/MSBradshaw/BOCC/main/WebPage/all_predicted_sign_clusters.2021.txt';
// an array of arrays where each sub-array is a member list of HPO terms and gene symbols, except index-0 is the cluster name
var clusters = [];

fetch(url)
  .then(function(response) {
    response.text().then(function(text) {
        var storedTextArr = text.split(/\r?\n/);
        for (var i = 0; i < storedTextArr.length; i++) {
            var tmp_members = storedTextArr[i].split(/\t/)
            clusters.push(tmp_members)
        }
    });
  });

function clean_and_split_input(text){
    var text_list = text.split(/\n/);
    // remove white spaces from HPO terms
    var clean_list = [];
    for (var i = 0; i < text_list.length; i++){
        item = text_list[i];
        if (typeof item == 'undefined'){
            continue;
        }
        item = item.replace(/\s/g, "");
        clean_list.push(item);
    }
    return clean_list;
}

function combine_and_bolden_cluster(gene,hpo,cluster){
    /*
    combine the cluster list into a string and place bold tags around the gene and HPO that matched in the cluster
    */
    var cluster_text = ""
    for(var i = 0; i < cluster.length; i++){
        if(i == 0){
            cluster_text += '<b>' + cluster[i] + '</b>';
            continue;
        }
        if( cluster[i] == gene || cluster[i] == hpo){
            cluster_text += ', <b>' + cluster[i] + '</b>';
        }else{
            cluster_text += ', ' + cluster[i];
        }
    }
    return cluster_text;
}

function search_clusters(hpos, genes, clusters){
    /*
    Input:
    Take a list of genes, hpos and clusters, display clusters containing gene & hpo co-occurances
    genes: list of genes. ['string1','string2'...]
    hpos: list of hpos. ['string1','string2'...]
    clusters: list clusters. [ ['gene1','hpo1','hpo2'...], ['gene2','gene3','hpo3'...], ... ]
    */
    var res_text = "";
    for (var i = 0; i < genes.length; i++) {
        gene = genes[i]
        for (var j = 0; j < hpos.length; j++) {
            hpo = hpos[j]
            for (var k = 0; k < clusters.length; k++) {
                cluster = clusters[k]
                if(cluster.includes(hpo, 1) & cluster.includes(gene, 1)){
                    console.log('Found')
                    console.log(hpo)
                    console.log(gene)
                    console.log(cluster)
                    res_text += '<div>' + combine_and_bolden_cluster(gene,hpo,cluster) + '</div></br>'
                }
            }
        }
    }
    document.getElementById("resultslist").innerHTML = res_text;
}

window.onload = function(){
    document.getElementById('search_button').onclick = function() {
        console.log('Clicked button!')
        // get the values stored in HPO section
        hpos = document.getElementById("hpolist").value
        var clean_hpos = clean_and_split_input(hpos);

        // get the values stored in the gene section
        genes = document.getElementById("genelist").value;
        var clean_genes = clean_and_split_input(genes);

        // search for co-occurances
        search_clusters(clean_hpos, clean_genes, clusters)

    };
};
