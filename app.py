from flask import Flask, request, jsonify, render_template
from Bio import Entrez


app = Flask(__name__)

# 设置Entrez的电子邮件（NCBI要求）
Entrez.email = "tony.liu.abm@gmail.com"


# 首页路由，显示表单
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/search_accession', methods=['POST'])
def search_accession():
    try:
        accession = request.form.get('accession')
        gene_id   = request.form.get('gene_id')

        print(accession, gene_id)
        if gene_id != '':
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
        elif accession != '':
            handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        data = handle.read()
        print(data)
        handle.close()
        return jsonify({'status': 'ok', 'data': data})
    except:
        return jsonify({'status': 'error', 'data': 'Failed to retrieve data'})
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5000,debug=True)
