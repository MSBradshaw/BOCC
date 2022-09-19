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

window.onload = function(){
    document.getElementById('search_button').onclick = function() {
        console.log('Clicked button!')
        // get the values stored in HPO section
        // get the values stored in the gene section
        // search for co-occurances
        for (var i = 0; i < clusters.length; i++) {
            console.log(clusters[i]);
            console.log('fucking fuck javascript')
        }
    };
};
