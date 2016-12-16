/**
 * Created by calarose on 28/11/16.
 */

window.onload = function() {

    /** Create the download link and download CSV file **/
    downloadCSV = function(csv, filename) {
        var csvFile = new Blob([csv], {
            type: 'text/csv'
        });

        var link = document.createElement('a');
        link.download = filename;
        link.href = window.URL.createObjectURL(csvFile);
        link.style.display = 'none';
        document.body.appendChild(link);
        link.click();
    };

    /** Generate the CSV file **/
    toCSV = function(filename) {
        var c, col, cols, file, i, j, len, len1, r, row, rows;
        rows = document.querySelectorAll('table tr');
        file = [];
        for (i = 0; i < rows.length; i++) {
            r = rows[i];
            row = [];
            cols = r.querySelectorAll('th, td');
            if (cols[0].innerText == ''){
                col = cols[0].innerText.split('\n').join(';');
                row.push(col);
                for(j = 1; j< cols.length; j++){
                    col = cols[j].innerText;
                    var colgroup = document.getElementById(col).getAttribute('colspan');
                    row.push(col);
                    for(var k = 1; k < colgroup; k++){
                        row.push('');
                    }
                }
            }
            else {
                for (j = 0; j < cols.length; j++) {
                    c = cols[j];
                    col = c.innerText;
                    col = col.split('\n').join('; ');
                    row.push(col);
                }
            }
            file.push(row.join(','));
        }
        downloadCSV(file.join('\r\n'), filename);
    };

    /** Generate the JSON file and download link, then download the file. **/
    downloadJSON = function (filename) {

        json = document.getElementsByTagName('p')[0].innerHTML;
        var jsonFile = new Blob([json], {
           type: 'application/json'
        });

        link = document.createElement('a');
        link.download = filename;
        link.href = window.URL.createObjectURL(jsonFile);
        link.style.display = 'none';
        document.body.appendChild(link);
        link.click();

    };


    if (document.getElementById('download-button-json') != null){
        document.querySelector('#download-button-json').addEventListener('click', function(e) {
            return downloadJSON('Results_Summary.json')
        });
    }
    if (document.getElementById('download-button-csv') != null){
        document.querySelector('#download-button-csv').addEventListener('click', function(e) {
            return toCSV('Results_Summary.csv');
        });
    }

};