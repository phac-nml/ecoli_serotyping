/**
 * Created by calarose on 28/11/16.
 */

window.onload = function() {

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

    toCSV = function(filename) {
        var c, col, cols, file, i, j, len, len1, r, row, rows;
        rows = document.querySelectorAll('table tr');
        file = [];
        for (i = 0; i < rows.length; i++) {
            r = rows[i];
            row = [];
            cols = r.querySelectorAll('th, td');
            for (j = 0; j < cols.length; j++) {
                c = cols[j];
                col = c.innerText;
                col = col.split('\n').join('; ');
                row.push(col);
            }
            file.push(row.join(','));
        }
        download(file.join('\r\n'), filename);
    };

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


    document.querySelector('#download-button-json').addEventListener('click', function(e) {
       return downloadJSON('Results_Summary.json')
    });
    document.querySelector('#download-button-csv').addEventListener('click', function(e) {
        return toCSV('Results_Summary.csv');
    });

};